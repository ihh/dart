#!/usr/bin/env perl -w

use Stockholm;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname <OPTIC multiple alignment file>    (or 'cat file | $progname')\n";
die $usage if @ARGV > 1;

my %stock;
while (<>) {
    next if /^group_id\s+gene\s+alignment\s*$/;
    chomp;
    my ($group_id, $gene, $alignment) = split;
    unless (exists $stock{$group_id}) {
	$stock{$group_id} = Stockholm->new();
	$stock{$group_id}->add_gf ("ID", $group_id);
    }
    $stock{$group_id}->add_row ($gene, $alignment);
}

for my $gene_id (sort { $a <=> $b } keys %stock) {
    print $stock{$gene_id}->to_string;
}
