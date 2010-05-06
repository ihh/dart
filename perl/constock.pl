#!/usr/bin/env perl -w

my $usage = "Usage: $0 [-tcs|-sps] [<Stockholm file(s)...>]\n\n";
$usage .= "Finds the \"consensus\" of a set of Stockholm-format multiple alignments,\n";
$usage .= "each of which is a global alignment of the same set of underlying sequences.\n";
$usage .= "Use -tcs for \"Total Column Score\", -sps for \"Sum-of-Pairs Score\" (SPS is default).\n";

my $sps = 1;
my @arg;
foreach my $arg (@ARGV) {
    if ($arg eq "-sps") { $sps = 1 }
    elsif ($arg eq "-tcs") { $sps = 0 }
    elsif ($arg =~ /^-/) { die $usage }
    else { push @arg, $arg }
}
@ARGV = @arg;

my $n_align = 0;
my (@seqname, %seq, @ungapped, %count, %rowindex, $rows);
while (<>) {
    if (/^\#/) {
	next;

    } elsif (/^\s*(\S+)\s+(\S+)\s*$/) {
	# alignment row
	my ($name, $seq) = ($1, $2);
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

    } elsif (/^\s*\/\/\s*$/) {
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
	for (my $col = 0; $col < $cols; ++$col) {
	    my @col = map (substr ($_, $col, 1), @row);  # get column as ASCII chars
	    my $n_ungapped = 0;
	    grep ($col[$_] eq "." || (++$n_ungapped && ++$pos[$_]), @rowindex);  # increment sequence coords & count non-gaps
	    next if $sps && $n_ungapped < 2;  # skip columns with zero or one ungapped chars
	    my $id = join ("", @col) . " @pos";  # unique text ID for column
	    $count{$id} = 0 unless defined $count{$id};
	    ++$count{$id};  # increment count for this column ID
	}

	# clear alignment info
	%seq = ();
	++$n_align;

    } elsif (/\S/) {  # back to main parse loop
	warn "Ignoring line: $_";
    }
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
warn "[Dynamic programming over ", @id+0, " column coords]\n";
my @score = map (0, @id);
my @max_score = @score;
my @traceback_col = map (undef, @id);
my $max_score = 0;
my $max_col;
for (my $dpcol = 0; $dpcol < $dpcols; ++$dpcol) {
    my ($col_text, $pos_array) = parse_column_id ($id[$dpcol]);
    my $score = 0;
    my $traceback_col;
    # find highest-socring previous compatible column
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
my @trace;
my $trace_col = $max_col;
push @trace, [map (length($_), @ungapped)];  # final column
while (defined $trace_col) {
    my ($col_text, $pos_array) = parse_column_id ($id[$trace_col]);
    push @trace, $pos_array;
    $trace_col = $traceback_col[$trace_col];
}
push @trace, [map (0, @ungapped)];  # initial column

# prepare output
my @output = map ("", @seqname);
my @pos = map (0, @seqname);
my @inc = @pos;
for (my $trace_index = @trace - 1; $trace_index >= 0; --$trace_index) {
    my @next_pos = @{$trace[$trace_index]};
    my $max_inc = 0;
    for (my $row = 0; $row < $rows; ++$row) {
	$inc[$row] = $next_pos[$row] - $pos[$row];
	$inc[$row] = 0 if $inc[$row] < 0;
	$max_inc = $inc[$row] if $max_inc < $inc[$row];
    }
    for (my $row = 0; $row < $rows; ++$row) {
	$output[$row] .= "." x ($max_inc - $inc[$row]) . substr ($ungapped[$row], $pos[$row], $inc[$row]);
    }
    @pos = @next_pos;
}

# print output
my $name_len = 0;
foreach my $name (@seqname) {
    $name_len = length($name) if length($name) > $name_len;
}
print "# STOCKHOLM 1.0\n";
for (my $row = 0; $row < $rows; ++$row) {
    print $seqname[$row], " " x ($name_len + 1 - length ($seqname[$row])), $output[$row], "\n";
}
print "//\n";

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
