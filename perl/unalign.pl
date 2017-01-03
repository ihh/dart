#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm;
use Newick;
use SequenceIterator qw(iterseq);

# command-line switches
my $randomize;

my $usage = "\n$0\n\nGiven a FASTA file, make a dummy Stockholm alignment(+tree),\nwherein all sequences are in fact unaligned.\n\nUsage: $0 [-rnd] <FASTA file>\n\nOptions:\n -rnd  Randomize input order (& hence tree)\n\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^(-h|-help|--help)$/) {
	die $usage;
    } elsif ($arg eq "-rnd") {
	$randomize = 1;
    } else {
	push @argv, $arg;
    }
}
unless (@argv) {
    @argv = ('-');
    warn "[waiting for alignments on standard input]\n";
}
die $usage unless @argv == 1;
my ($fasta) = @argv;

# read FASTA
my (@name, @seq, @length);
iterseq ($fasta,
	 sub {
	     my ($name, $seq) = @_;
	     $name =~ s/ .*//;
	     push @name, $name;
	     push @seq, $seq;
	     push @length, length $seq;
	 });

# randomize?
if ($randomize) {
    # Fisher-Yates shuffle
    my @order = (0 .. @name - 1);
    for my $i (0 .. @name - 2) {
	my $j = $i + int (rand (@name - $i));
	if ($i != $j) {
	    @order[$i,$j] = @order[$j,$i];
	}
    }
    @name = @name[@order];
    @seq = @seq[@order];
    @length = @length[@order];
}

# calculate indentation offsets
my $offset = 0;
my @offset = (0, map ($offset += $_, @length));
pop @offset;

# create alignment
my $stock = Stockholm->new;
for my $i (0..@name-1) {
    $stock->add_row ($name[$i], ("." x $offset[$i]) . $seq[$i] . ("." x ($offset - $offset[$i] - $length[$i])));
}

my @parent = map (-1, @name);
for (my $i = 0; $i < @parent - 1; $i += 2) {
    my $next_node = @parent;
    $parent[$i] = $parent[$i+1] = $next_node;
    push @parent, -1;
    push @name, undef;
}

my @rev_name = reverse @name;
my @rev_parent = map ($_ < 0 ? $_ : (@parent - 1 - $_), reverse @parent);

my $tree = Newick->new;
@{$tree->parent} = @rev_parent;
@{$tree->node_name} = @rev_name;
@{$tree->branch_length} = (undef, map (1, 1..@rev_name));

$stock->gf_NH ([$tree->to_string]);

print $stock->to_string;
