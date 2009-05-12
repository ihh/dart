#!/usr/bin/perl -w

use Stockholm;
use Newick;
use SequenceIterator qw(iterseq);

unless (@ARGV) {
    @ARGV = qw(-);
    warn "[waiting for FASTA file on standard input]\n";
}
die "\n$0\n\nGiven a FASTA file, make a dummy Stockholm alignment,\nwherein all sequences are in fact unaligned.\n\nUsage: $0 <FASTA file>\n\n"
    if @ARGV != 1 || grep (/^(-h|-help|--help)$/, @ARGV);
my ($fasta) = @ARGV;

my (@name, @seq, @offset);
my $offset = 0;
iterseq ($fasta,
	 sub {
	     my ($name, $seq) = @_;
	     $name =~ s/ .*//;
	     push @name, $name;
	     push @seq, $seq;
	     push @offset, $offset;
	     $offset += length $seq;
	 });

my $stock = Stockholm->new;
for my $i (0..@name-1) {
    $stock->add_row ($name[$i], ("." x $offset[$i]) . $seq[$i] . ("." x ($offset - $offset[$i] - length($seq[$i]))));
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
