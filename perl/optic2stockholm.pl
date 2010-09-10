#!/usr/bin/env perl -w

use Stockholm;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $tree = "nj";
my $seq = "gene";
my $type = "aa";

my $usage = "Usage:\n $progname [-nj|-ds] [-gene|-transcript] [-aa|-na] <OPTIC clade dir>\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-(nj|ds)$/) {
	$tree = $1;
    } elsif ($arg =~ /^-(gene|transcript)$/) {
	$seq = $1;
    } elsif ($arg =~ /^-(aa|na)$/) {
	$type = $1;
    } elsif ($arg =~ /^-/) {
	die $usage;
    } else {
	push @argv, $arg;
    }
}

die $usage if @argv != 1;
my ($opticRoot) = @argv;

my %stock;
my $alignments = 0;
local *MULTI;
open MULTI, "gzip -cd $opticRoot/multiple_alignments/${seq}s_${type}.gz |" or die $!;
while (<MULTI>) {
    next if /^group_id\s+(gene|transcript)\s+alignment\s*$/;
    chomp;
    my ($group_id, $name, $alignment) = split;
    unless (exists $stock{$group_id}) {
	warn "Creating alignment #$group_id\n" if (++$alignments) % 100 == 0;
	$stock{$group_id} = Stockholm->new();
	$stock{$group_id}->add_gf ("ID", $group_id);
    }
    $stock{$group_id}->add_row ($name, $alignment);
}
close MULTI or die $!;

my $trees = 0;
my $group_id;
local *TREE;
open TREE, "gzip -cd $opticRoot/orthologs/${tree}_trees.gz |" or die $!;
while (<TREE>) {
    if (/^>(\S+)/) {
	$group_id = $1;
    } elsif (/\S/ && defined($group_id) && exists $stock{$group_id}) {
	chomp;
	$stock{$group_id}->add_gf ("NH", $_);
	warn "Attached tree #$group_id\n" if (++$trees) % 100 == 0;
    }
}
close TREE or die $!;

for my $group_id (sort { $a <=> $b } keys %stock) {
    print $stock{$group_id}->to_string;
}
