#!/usr/bin/perl

use Getopt::Long;

use Newick;
use Stockholm;
use Stockholm::Database;

my $tag = "NH";
my $topo = 0;

my $usage = "";
$usage .= "$0 -- extract a Newick tree annotation from a Stockholm alignment\n";
$usage .= "\n";
$usage .= "Usage: $0 [-gf <tag>] [-topology] [filename(s)]\n";
$usage .= "\n";
$usage .= "The 'tag' is the Stockholm '#=GF XX' tag that the tree is stored under.\n";
$usage .= "By default it is '$tag'.\n";
$usage .= "\n";

GetOptions ("tag=s" => \$tag, "topology" => \$topo) or die $usage;

unless (@ARGV) {
    @ARGV = qw(-);
    warn "[waiting for alignments on standard input]\n";
}

while (@ARGV) {
    my $db = Stockholm::Database->from_file (shift);
    for my $stock (@$db) {
	my $newick = $stock->newick ($tag);
	if (defined $newick) {
	    if ($topo) { grep ($_ = undef, @{$newick->branch_length}) }
	    print $newick->to_string, "\n";
	}
    }
}
