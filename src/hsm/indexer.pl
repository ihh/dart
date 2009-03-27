#!/usr/local/bin/perl -w

my $pwd = `pwd`;
chomp $pwd;

my $pfam_dir = "Pfam";
my $index_filename = "Index";
my $align_dir = "$pfam_dir/Alignments";
my $align_filename = "Seed.mul";

my $tree_file = "Tree";

local *INDEX;
open INDEX, "<$pfam_dir/$index_filename" or die $!;
while (<INDEX>) {
    my ($dir, $name) = split;
    $dir = "$align_dir/$dir";
    print "$name $dir/$tree_file $dir/$align_filename\n";
}
close INDEX;
