#!/usr/bin/perl -w

use Stockholm;

@ARGV = qw(-) unless @ARGV;
die "Usage: $0 <Stockholm_file>" unless @ARGV == 1;

my ($filename) = @ARGV;
my $stock = Stockholm->from_file ($filename);

my $cols = $stock->columns;
my @gff = @{$stock->gf_GFF};

my @rmean = map (0, 1..$cols);
my $maxState = 0;
for my $gff (@gff) {
    my @f = split /\t/, $gff;
    if ($f[4] == $cols) {
	if ($f[8] =~ /lgPost=([\-\d\.\+\-eE]+)/) {
	    my $lgPost = $1;
	    if ($f[2] =~ /^S(\d+)$/) {
		my $state = $1;
		$maxState = $state if $state > $maxState;
		$rmean[$f[3]] += $state * 2**$lgPost;
	    }
	}
    }
}

print "track type=wiggle_0\n";
print "fixedStep chrom=XrateAlignment start=1 step=1\n";
for my $pos (0..$#rmean) {
    print $rmean[$pos]/$maxState, "\n";
}

