#!/usr/local/bin/perl -w

my $pwd = `pwd`;
chomp $pwd;

my $pfam_dir = "Pfam";
my $index_filename = "Index";
my $align_dir = "$pfam_dir/Alignments";
my $align_filename = "Seed.mul";

my $logfile = "WEIGHPFAM.LOG";

my $smdist_prog = "/users/ihh/bin/smdist -log 6 -logfile SMDIST.LOG";
my $weighbor_prog = "/users/ihh/bin/weighbor -v";
my $crapfilter_prog = "$pwd/crapfilter.pl";

my $distance_matrix_file = "Distances";
my $tree_file = "Tree";

local *ALIGN_DIR;
opendir ALIGN_DIR, $align_dir or die $!;
my @align = grep -d, map "$align_dir/$_", grep !/^\./, readdir ALIGN_DIR;
closedir ALIGN_DIR;

open LOGFILE, ">$logfile" or die $!;
foreach my $align (@align)
{
    my $command = "cd $align; $smdist_prog $align_filename |$crapfilter_prog >$distance_matrix_file; $weighbor_prog -i $distance_matrix_file -o $tree_file";
    print LOGFILE $command, "\n";
    system $command;
}
close LOGFILE or die $!;
