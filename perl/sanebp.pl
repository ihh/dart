#!/usr/bin/env perl -w

use FileHandle;

my $usage = "Usage: $0 <Stockholm database>\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } else {
	push @argv, $arg;
    }
}
die $usage unless @argv == 1;
my ($stockfile) = @argv;

warn "Reading alignments from $stockfile\n";
my $in = FileHandle->new;
$in->open ("<$stockfile") or die "Couldn't open $stockfile: $!";

local $_;
while ($_ = $in->getline()) {
    next unless /\S/;
    if (/^\s*(\S+)\s+(\S+)\s*$/) {
	my ($name, $seq) = ($1, $2);
	$stock = addTag ($stock, $name, "", $seq);
    } elsif (/^\s*\#=GR\s*(\S+)\s*(\S+)\s*(\S+)/) {
	my ($name, $tag, $data) = ($1, $2, $3);
	$stock = addTag ($stock, $name, $tag, $data);
    } elsif (/\s*\/\/\s*$/) {
	sanitise ($stock);
	$stock = undef;
    } else {
	$stock = addTag ($stock, "", "", $_);
    }
}
sanitise ($stock) if defined $stock;
$in->close();

sub is_lchar {
    my ($s) = @_;
    return $s eq '<' || $s eq '(' || $s eq '[' || $s eq '{';
}

sub is_rchar {
    my ($s) = @_;
    return $s eq '>' || $s eq ')' || $s eq ']' || $s eq '}';
}

sub is_gap {
    my ($s) = @_;
    return $s eq '-' || $s eq '.' || $s eq '_';
}

sub sanitise {
    my ($stock) = @_;
    my @seqs = grep (length>0, keys %$stock);
    return if @seqs == 0;
    my $firstseqname = $seqs[0];
    my $firstseq = $stock->{$firstseqname}->{""};
    my $cols = length ($firstseq);
    my @row = map (\$stock->{$_}->{""}, @seqs);
    my @ss = map (\$stock->{$_}->{"SS"}, @seqs);
    for (my $row = 0; $row < @row; ++$row) {
	my $name = $seqs[$row];
	my $r = $row[$row];
	my $ss = $ss[$row];
	if (length($$r) != $cols) { die "Sequence $name has length ", length($$r), "; expected $cols" }
	if (length($$ss) != $cols) { die "SS line of sequence $name has length ", length($$ss), "; expected $cols" }
	my @lcol;
	my $removed = 0;
	for (my $col = 0; $col < $cols; ++$col) {
	    my $s = substr ($$ss, $col, 1);
	    if (is_lchar($s)) { push @lcol, $col }
	    elsif (is_rchar($s)) {
		if (@lcol == 0) { die "Sequence $name has bad fold string: too many >'s in $$ss" }
		my $lcol = pop @lcol;
		# check that all rows either have gaps or secondary structure at both positions
		for (my $i = 0; $i < @row; ++$i) {
		    if ($i != $row) {
			my $lc = substr (${$row[$i]}, $lcol, 1);
			my $rc = substr (${$row[$i]}, $col, 1);
			my $ls = substr (${$ss[$i]}, $lcol, 1);
			my $rs = substr (${$ss[$i]}, $col, 1);
			my $both_gap = is_gap($lc) && is_gap($rc) && is_gap($ls) && is_gap($rs);
			my $both_ss = !is_gap($lc) && !is_gap($rc) && is_lchar($ls) && is_rchar($rs);
#			warn "[checking basepair (", $lcol+1, ",", $col+1, ") in sequence $name: row=$row i=$i lc=$lc rc=$rc ls=$ls rs=$rs]\n";
			# if not gap or SS, then remove this basepair
			unless ($both_gap || $both_ss) {
			    substr ($$ss, $lcol, 1) = '.';
			    substr ($$ss, $col, 1) = '.';
			    ++$removed;
			    last;
			}
		    }
		}
	    }
	}
	if (@lcol != 0) { die "Sequence $name has bad fold string: too many <'s in $$ss" }
	if ($removed) { warn "[removed $removed basepairs from sequence $name]\n" }
    }
    # print out the alignment (without much care)
    while (my ($name, $taghash) = each %$stock) {
	while (my ($tag, $val) = each %$taghash) {
	    if (length $tag) { print "#=GR $name $tag $val\n" }
	    elsif (length $name) { print "$name $val\n" }
	    else { print $val }
	}
    }
    print "//\n";
}

sub addTag {
    my ($stock, $name, $tag, $data) = @_;
    if (!defined ($stock)) {
	$stock = {};
    }
    if (!exists $stock->{$name}) {
	$stock->{$name} = {};
    }
    if (!exists $stock->{$name}->{$tag}) {
	$stock->{$name}->{$tag} = "";
    }
    $stock->{$name}->{$tag} .= $data;
    return $stock;
}
