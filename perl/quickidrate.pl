#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

my $usage = "\nUsage: $0 <total length> <number of sequences> [gap frequency file]\n\nEstimates insertion and deletion rates.\nExpects three-column input format, one line per pairwise alignment:\n\nColumn #1 is the 'substitution distance' between the two sequences\nColumn #1 is the total number of match columns in the alignment\nColumn #3 is the number of match columns which are directly followed by match columns\n\n";

die $usage if @ARGV < 2 || @ARGV > 4 || grep /^-/, @ARGV;

my $total = shift;
my $seqs = shift;

my $eqm_extend_prob = $total / ($total + $seqs);

my $dwell_time = 0;
my $events = 0;
while (<>) {
    if (/\d/) {
	my ($dist, $m, $mm) = split;
	$dwell_time += $dist * $m;
	$events += $m - $mm;
    }
}

my $delete_rate = $events / $dwell_time;
my $insert_rate = $delete_rate * $eqm_extend_prob;

print "Insert rate $insert_rate\nDelete rate $delete_rate\n";

warn "Umm, total length < number of sequences. Did you get them the right way round?\n" if $total < $seqs;


