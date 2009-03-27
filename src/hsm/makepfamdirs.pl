#!/usr/local/bin/perl -w

my $pfam_dir = "Pfam";
my $index_filename = "Index";
my $align_dir = "$pfam_dir/Alignments";
my $align_filename = "Seed.mul";

mkdir $pfam_dir, 0777 unless -e $pfam_dir;
mkdir $align_dir, 0777 unless -e $align_dir;

my @index;
my ($ac, $id, @align);
while (<>) {
    if (/^\/\//) {
	my $align_subdir = "$align_dir/$ac";
	mkdir $align_subdir, 0777 unless -e $align_subdir;
	local *ALIGN;
	open ALIGN, ">$align_subdir/$align_filename" or die $!;
	print ALIGN @align;
	close ALIGN or die $!;
	push @index, "$ac $id\n";
	warn "$ac $id\n";
	($ac, $id, @align) = ();
    } elsif (/^\#=GF\s+ID\s+(.+)/) { $id = $1; $id =~ s/\s/_/g }
    elsif (/^\#=GF\s+AC\s+(\S+)/) { $ac = $1 }
    elsif (!/^\#/) { push @align, $_ }
}
local *INDEX;
open INDEX, ">$pfam_dir/$index_filename" or die $!;
print INDEX @index;
close INDEX or die $!;
