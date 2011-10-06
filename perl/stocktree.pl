#!/usr/bin/perl

use Getopt::Long;
use Newick;
use Stockholm;

my $tag = "NH";

my $usage = "";
$usage .= "$0 -- extract a Newick tree annotation from a Stockholm alignment\n";
$usage .= "\n";
$usage .= "Usage: $0 [-gf <tag>] [filename(s)]\n";
$usage .= "\n";
$usage .= "The 'tag' is the Stockholm '#=GF XX' tag that the tree is stored under.\n";
$usage .= "By default it is '$tag'.\n";
$usage .= "\n";

GetOptions ("tag=s" => \$tag) or die $usage;

unless (@ARGV) {
    @ARGV = qw(-);
    warn "[waiting for alignments on standard input]\n";
}

while (@ARGV) {
    my $stock = Stockholm->from_file (shift);
    print map ("$_\n", @{$stock->gf->{$tag}});
}
