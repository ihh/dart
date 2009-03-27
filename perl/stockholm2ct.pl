#!/usr/bin/perl -w

my $PS = "PS";
my $SS = "SS";
my $SC = "SC";
my $PS_CONS = "PS_cons";
my $SS_CONS = "SS_cons";

my $ct_prefix = "STRUCT_";
my $ct_number_prefix = "_ALIGN_";
my $ct_suffix = ".ct";

my %rchar = ('('=>')','<'=>'>','['=>']');
my %lchar = map (($rchar{$_}=>$_), keys(%rchar));
my %gapchar = map (($_=>1), '-', '.', '_', ',');

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $pad = 0;

my $usage = "Usage: $progname <Stockholm file>\n";
$usage .=   "             [-h] print this help message\n";
$usage .=   "    [-p <prefix>] change prefix for CT files (default is '$ct_prefix')\n";
# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-p") {
	defined ($ct_prefix = shift) or die $usage;
    } else {
	push @argv, $arg;
    }
}
die $usage unless @argv == 1;

my ($stockfile) = @argv;

local *STOCK;
open STOCK, "<$stockfile" or die "Couldn't open '$stockfile': $!";

my $ct_number = 0;
while (1) {
    my (@seqname, %seq, %ss, $ps_cons, $ss_cons, $score);
    my $found_separator = 0;
    while (<STOCK>) {
	if (/^\s*\#=GR\s+(\S+)\s+$SS\s+(\S+)\s*$/) { $ss{$1} .= $2 }  # by-seq secondary structure
	elsif (/^\s*\#=GC\s+$PS_CONS\s+(\S+)\s*$/) { $ps_cons .= $1 }  # consensus primary sequence
	elsif (/^\s*\#=GC\s+$SS_CONS\s+(\S+)\s*$/) { $ss_cons .= $1; }  # consensus secondary structure
	elsif (/^\s*\#=GF\s+$SC\s+(\S+)\s*$/) { $score = $1 }  # score
	elsif (/^\s*\#/) { }  # unrecognised line starting with '#'; do nothing
	elsif (/^\s*\/\//) { $found_separator = 1; last }  # alignment separator
	elsif (/^\s*(\S+)\s*(\S+)\s*$/) { $seq{$1} .= $2; push @seqname, $1 }
    }
    if (@seqname) {
	warn "Creating CT files for alignment #", ++$ct_number, "\n";
    }
    # loop through sequence names
    foreach my $encoded_seqname (@seqname) {
	# get alignment row data & sequence length
	my $row = $seq{$encoded_seqname};
	my $ss = exists($ss{$encoded_seqname}) ? $ss{$encoded_seqname} : $ss_cons;
	my $cols = length $row;
	# sequence name may encode co-ords of local alignment
	my ($seqname, $start, $end);
	if ($encoded_seqname =~ /^(\S+)\/(\d+)\-(\d+)$/) {
	    ($seqname, $start, $end) = ($1, $2, $3);
	} else {
	    ($seqname, $start, $end) = ($encoded_seqname, 1, undef);
	}
	# get basepair co-ords
	die "Sequence $encoded_seqname and its fold string have different lengths ($cols,", length($ss), ")" unless length($ss) == $cols;
	my @ungapped_pos = map (-1, 1..$cols);
#	my $ungapped = 'N' x ($start-1);
	my $ungapped = "";
	my @pair_pos = map (-1, 2..$start);
	my @char_stack;
	my @lpos_stack;
	for (my $i = 0; $i < $cols; ++$i) { # loop through columns in alignment, left to right
	    my $c = substr ($ss, $i, 1);  # get next char in fold string
	    my $seqchar = substr ($row, $i, 1);  # get next char in alignment row
	    my $pos;
	    if (exists $gapchar{$seqchar}) {  # check if this sequence doesn't have a gap in this column of the alignment
		$pos = -1;
	    } else {
		$pos = length $ungapped;  # keep track of the index into the ungapped sequence
		$ungapped .= $seqchar;
		push @pair_pos, -1;
	    }
	    $ungapped_pos[$i] = $pos;
	    if (exists $rchar{$c}) {  # if char is a '<' or equivalent, push char & index onto respective stacks
		#warn "c=$c, pushing $pos";
		push @lpos_stack, $pos;
		push @char_stack, $c;
	    } elsif (exists $lchar{$c}) {  # if char is a '>' or equivalent then...
		if (@char_stack==0 || pop(@char_stack) ne $lchar{$c}) {  # ...pop a char off the char stack, check it matches
		    die "Position $i of $seqname: unmatched $c\n";
		}
		my $lpos = pop @lpos_stack;  # pop index off the index stack
		#warn "c=$c, pos=$pos, popping $lpos";
		if ($lpos >= 0 && $pos >= 0) {  # check that neither of the basepaired residues is actually a gap
		    # pair up the two residues
		    $pair_pos[$lpos] = $pos;
		    $pair_pos[$pos] = $lpos;
		}
	    }
	}
	my $seqlen = length $ungapped;  # calculate length of sequence
	my $coordlen = $start > $end ? ($start + 1 - $end) : ($end + 1 - $start);
	if ($coordlen != $seqlen) {
		warn "Sequence co-ords in '$encoded_seqname' inconsistent with apparent sequence length of $seqlen\n";
	    }
	# create CT filename
	my $filename = $seqname;
	$filename =~ s/\//_/g;
	$filename = $ct_prefix . $filename . $ct_number_prefix . $ct_number . $ct_suffix;
	# write CT file
	warn " (creating file '$filename')\n";
	local *CTFILE;
	open CTFILE, ">$filename" or die "Couldn't open '$filename': $!";
	print CTFILE "$seqlen ENERGY = 0    $seqname    $seqlen bp RNA\n";
	for (my $pos = 0; $pos < $seqlen; ++$pos) {
	    printf CTFILE "%3d %s %7d %4d %4d %4d\n", $pos+1, substr($ungapped,$pos,1), $pos, $pos+2, $pair_pos[$pos]+1, $pos+1;
	}
	close CTFILE or die "Couldn't close '$filename': $!";
    }
    # if no alignment separator, then quit the loop
    last unless $found_separator;
}
close STOCK;
