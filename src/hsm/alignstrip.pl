#!/usr/local/bin/perl -w

my $pwd = `pwd`;
chomp $pwd;

my $pfam_dir = "Pfam";
my $index_filename = "Index";
my $align_dir = "$pfam_dir/Alignments";
my $align_file = "Seed.mul";

my $crapfilter_prog = "$pwd/crapfilter.pl";

local *ALIGN_DIR;
opendir ALIGN_DIR, $align_dir or die $!;
my @align = grep -d, map "$align_dir/$_", grep !/^\./, readdir ALIGN_DIR;
closedir ALIGN_DIR;

foreach my $align (@align)
{
    my $command = "cd $align; mv $align_file $align_file.OLD; cat $align_file.OLD |$crapfilter_prog >$align_file; rm $align_file.OLD";
    warn $command, "\n";
    system $command;
}
