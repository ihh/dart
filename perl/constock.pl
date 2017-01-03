#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm;
use Stockholm::Database;
use Newick;

my $postprob_tag = 'PP';

my $usage = "Usage: $0 [-tcs|-sps] [-pp|-aa] [<Stockholm file(s)...>]\n\n";
$usage .= "Finds the \"consensus\" of a set of Stockholm-format multiple alignments,\n";
$usage .= "each of which is a global alignment of the same set of underlying sequences.\n";
$usage .= "Use -tcs for \"Total Column Score\", -sps for \"Sum-of-Pairs Score\" (SPS is default).\n";
$usage .= "Use -pp for posterior-probability coloring in output, -aa for amino-acid coloring.\n";

my $sps = 1;
my $pipe_to_cat = 0;
my @arg;
foreach my $arg (@ARGV) {
    if ($arg eq "-sps") { $sps = 1 }
    elsif ($arg eq "-tcs") { $sps = 0 }
    elsif ($arg eq "-aa") { $color = "AMINO" }
    elsif ($arg eq "-pp") { $color = \&pp_color_scheme }
    elsif ($arg eq "-less") { $pipe_to_cat = 0 }
    elsif ($arg =~ /^-/) { die $usage }
    else { push @arg, $arg }
}
@ARGV = @arg;

unless (@ARGV) { @ARGV = ("-"); warn "[waiting for Stockholm file(s) on standard input]\n" }
my @db;
for my $filename (@ARGV) {
    push @db, @{Stockholm::Database->from_file($filename)};
}

my $n_align = 0;
my (@seqname, %seq, @ungapped, %count, %rowindex, $rows);
for my $stock (@db) {
    warn "...processing alignment #", $n_align+1, " of ", @db+0, "\n";
    my $tree;
    if (defined $stock->gf->{'NH'}) {
	$tree = Newick->from_string (@{$stock->gf_NH});
    }
    while (my ($name, $seq) = each %{$stock->seqdata}) {
	if (defined $tree) {
	    # drop internal nodes
	    my @node = $tree->find_nodes ($name);
	    next if @node == 1 && $tree->children($node[0]) > 0;
	}
	if (!exists $seq{$name}) {
	    if ($n_align == 0) {
		$rowindex{$name} = @seqname;
		push @seqname, $name;
	    } else {
		die "Previously unseen sequence in alignment ", $n_align+1, ": $name\n" unless exists $rowindex{$name};
	    }
	    $seq{$name} = "";
	}
	$seq = lc $seq;
	$seq =~ s/\-/\./g;  # convert "-" gap characters into "."
	$seq{$name} .= $seq;
    }

    # end of alignment: check consistency
    for (my $row = 0; $row < @seqname; ++$row) {
	my $name = $seqname[$row];
	die "Missing sequence in alignment ", $n_align+1, ": $name\n" unless exists $seq{$name};
	my $ungapped = $seq{$name};
	$ungapped =~ s/\.//g;
	if ($n_align == 0) {
	    $ungapped[$row] = $ungapped;
	} elsif ($ungapped ne $ungapped[$row]) {
	    die "Sequence '$name' in alignment ", $n_align+1,
	    " does not match earlier alignments\nNew: $ungapped\nOld: $ungapped[$row]\n";
	}
    }

    # count rows & columns, and (again) check consistency
    $rows = @seqname if $n_align == 0;
    my ($refname, $cols);
    while (my ($name, $seq) = each %seq) {
	my $l = length $seq;
	if (!defined $cols) {
	    $cols = $l;
	    $refname = $name;
	} elsif ($l != $cols) {
	    die "In alignment ", $n_align+1, ": row $name has length $l, while row $refname has length $cols\n";
	}
    }

    # get coords of each column
    my @rowindex = 0..$rows-1;  # useful array of row indices
    my @row = map ($seq{$_}, @seqname);  # sort alignment rows
    my @pos = map (0, @seqname);  # initialise sequence coords
    my $pp = $stock->gc->{$postprob_tag};
    for (my $col = 0; $col < $cols; ++$col) {
	my @col = map (substr ($_, $col, 1), @row);  # get column as ASCII chars
	my $n_ungapped = 0;
	grep ($col[$_] eq "." || (++$n_ungapped && ++$pos[$_]), @rowindex);  # increment sequence coords & count non-gaps
	next if $sps && $n_ungapped < 2;  # skip columns with zero or one ungapped chars
	my $id = join ("", @col) . " @pos";  # unique text ID for column
	my $weight = defined($pp) ? ((1 + substr($pp,$col,1)) / 10) : 1;
	$count{$id} = 0 unless defined $count{$id};
	$count{$id} += $weight;  # increment count for this column ID
    }

    # clear alignment info
    %seq = ();
    ++$n_align;
}

# if no alignments, bail here
exit unless @seqname;

# get rank of each column ID
my @id = keys %count;
my @rank;
foreach my $id (@id) {
    my ($col_text, $pos_array) = parse_column_id ($id);
    my $rank = 0;
    foreach my $pos (@$pos_array) { $rank += $pos }
    push @rank, $rank;
}

# sort IDs by rank
my $dpcols = @id;
@id = @id [sort { $rank[$a] <=> $rank[$b] } 0..$dpcols-1];
my @count = map ($count{$_}, @id);
@rank = %count = ();  # clear unused stuff

# main DP loop
warn "[Dynamic programming over $dpcols column coords]\n";
my @score = map (0, @id);
my @max_score = @score;
my @traceback_col = map (undef, @id);
my $max_score = 0;
my $max_col;
my $report_interval = int ($dpcols / @db) + 1;
for (my $dpcol = 0; $dpcol < $dpcols; ++$dpcol) {
    if ($dpcol % $report_interval == 0) {
	warn "...processing column ", $dpcol+1, " of ", $dpcols, "\n";
    }
    my ($col_text, $pos_array) = parse_column_id ($id[$dpcol]);
    my $score = 0;
    my $traceback_col;
    # find highest-scoring previous compatible column
    for (my $prev_dpcol = $dpcol - 1; $prev_dpcol >= 0 && $score < $max_score[$prev_dpcol]; --$prev_dpcol) {
	my $prev_score = $score[$prev_dpcol];
	if ($prev_score > $score && compatible ($col_text, $pos_array, parse_column_id ($id[$prev_dpcol]))) {
	    $score = $prev_score;
	    $traceback_col = $prev_dpcol;
	}
    }

    # add score for this column, equal to count*weight, where weight depends on whether "-sps" or "-tcs" was selected
    my $weight;
    if ($sps) {
	# SPS: weight is equal to count*N*(N-1), where N is the number of non-gap characters
	my $non_gap_chars = 0;
	for (my $row = 0; $row < $rows; ++$row) {
	    if (substr ($col_text, $row, 1) ne ".") {
		++$non_gap_chars;
	    }
	}
	$weight = $non_gap_chars * ($non_gap_chars - 1);
    } else {
	# TCS: weight is equal to 1
	$weight = 1;
    }
    $score += $count[$dpcol] * $weight;

    # update DP matrix
    if ($score > $max_score) {
	$max_score = $score;
	$max_col = $dpcol;
    }
    $score[$dpcol] = $score;
    $max_score[$dpcol] = $max_score;
    $traceback_col[$dpcol] = $traceback_col;
}

# traceback
my (@trace, @trace_score);
my $trace_col = $max_col;
push @trace, [map (length($_), @ungapped)];  # final column
push @trace_score, 1;
while (defined $trace_col) {
    my ($col_text, $pos_array) = parse_column_id ($id[$trace_col]);
    push @trace, $pos_array;
    push @trace_score, $count[$trace_col] / $n_align;
    $trace_col = $traceback_col[$trace_col];
}
push @trace, [map (0, @ungapped)];  # initial column
push @trace_score, 1;

# prepare output
my @output = map ("", @seqname);
my @pos = map (0, @seqname);
my @inc = @pos;
my $score_line;
for (my $trace_index = @trace - 1; $trace_index >= 0; --$trace_index) {
    my @next_pos = @{$trace[$trace_index]};
    my $max_inc = 0;
    for (my $row = 0; $row < $rows; ++$row) {
	$inc[$row] = $next_pos[$row] - $pos[$row];
	if ($inc[$row] < 0) { $inc[$row] = 0; $next_pos[$row] = $pos[$row] }
	$max_inc = $inc[$row] if $max_inc < $inc[$row];
    }
    for (my $row = 0; $row < $rows; ++$row) {
	$output[$row] .= "." x ($max_inc - $inc[$row]) . substr ($ungapped[$row], $pos[$row], $inc[$row]);
    }
    $score_line .= int ($trace_score[$trace_index] * 9) x $max_inc;
    @pos = @next_pos;
}

# print output
my $stock = Stockholm->new;
for (my $row = 0; $row < $rows; ++$row) {
    $stock->add_row ($seqname[$row], $output[$row]);
}
$stock->gc->{$postprob_tag} = $score_line;

my $pager = $pipe_to_cat ? "cat" : "less -r";
local *LESS;
open (LESS, "|$pager") or open (LESS, ">-");

my @score_col = ([0,4], [0,2], [0,1], [0,6], [0,3], [0,7], [4,7], [2,7], [3,7], [6,7]);
print LESS $stock->to_string ("MAXCOLS" => "SCREEN", defined($color) ? ("COLOR" => $color) : ());

close LESS;

sub pp_color_scheme {
    my ($residue, $seqname, $pos) = @_;
    my $trace_sc = substr ($score_line, $pos, 1);
    my @sc = @{$score_col[$trace_sc]};
    return @{$score_col[$trace_sc]};
}

# functions

sub parse_column_id {
    my ($id) = @_;
    my ($col_text, @pos) = split /\s/, $id;
    return ($col_text, \@pos);
}

sub compatible {
    my ($col_text, $pos_array, $prev_col_text, $prev_pos_array) = @_;
    my $rows = @$pos_array;
    for (my $row = 0; $row < $rows; ++$row) {
	if ($pos_array->[$row] <= $prev_pos_array->[$row]
	    && substr ($col_text, $row, 1) ne "."
	    && substr ($prev_col_text, $row, 1) ne ".") {
	    return 0;
	}
    }
    return 1;
}
