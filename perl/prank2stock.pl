#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Newick;
use Stockholm;

# usage message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "$progname: convert PRANK output to Stockholm-with-Newick\n";
$usage .=   "\n";
$usage .=   "Usage: $progname <PRANK output file> <original Newick file>\n";
$usage .=   "              [-h]  print this message\n";
$usage .=   "\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } else {
	push @argv, $arg;
    }
}
die $usage unless @argv == 2;
my ($prankfile, $newickfile) = @argv;


# read PRANK output
local *PRANK;
open PRANK, "<$prankfile" or die "$prankfile : $!";

my $treeline = <PRANK>;
$treeline =~ s/^#\s+(\S+).*$/$1;/;
my $newtree = Newick->from_string ($treeline);

my $stock = Stockholm->new();
my ($seqname);
while (<PRANK>) {
    if (/^>(\S+)/) {
	$seqname = $1;
    } elsif (/(\S+)/) {
	chomp;
	$stock->add_row ($seqname, $1);
    }
}

close PRANK;


# read original tree
my $oldtree = Newick->from_file ($newickfile);


# check for simple match between trees
if ($oldtree->nodes != $newtree->nodes) { die "number of nodes doesn't match" }
my @oldleaves = $oldtree->leaves;
my @newleaves = $newtree->leaves;
if (@oldleaves != @newleaves) { die "number of leaves doesn't match" }
if (grep ($oldleaves[$_] != $newleaves[$_], 0..$#oldleaves)) { die "leaf indices don't match" }
if (grep ($oldtree->node_name->[$_] ne $newtree->node_name->[$_], @oldleaves)) { die "leaf names don't match" }
if (grep ($oldtree->parent->[$_] != $newtree->parent->[$_], 0..$newtree->nodes-1)) { die "parent indices don't match" }

# create new->old seqname map
my %new2old;
for my $n (0..$newtree->nodes-1) {
    my $oldname = $oldtree->node_name->[$n];
    my $newname = $newtree->node_name->[$n];
    $new2old{$newname} = $oldname;
}

# apply old->new seqname map
@{$stock->seqname} = map ($new2old{$_}, @{$stock->seqname});
%{$stock->seqdata} = map (($new2old{$_} => $stock->seqdata->{$_}), keys %{$stock->seqdata});

# output
$stock->add_gf ("NH", $oldtree->to_string);
print $stock->to_string ($stock->columns + $stock->maxNameLen + 1);
