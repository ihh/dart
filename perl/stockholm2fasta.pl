#!/usr/bin/env perl -w

my $columns = 50;
my $gapped = 0;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname [<Stockholm file(s)>]\n";
$usage .=   "             [-h] print this help message\n";
$usage .=   "             [-g] write gapped FASTA output\n";
$usage .=   "             [-s] sort sequences by name\n";
$usage .=   "      [-c <cols>] number of columns for FASTA output (default is $columns)\n";
# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-g") {
	$gapped = 1;
    } elsif ($arg eq "-s"){
	$sorted = 1;
    } elsif ($arg eq "-c") {
	defined ($columns = shift) or die $usage;
    } else {
	push @argv, $arg;
    }
}
@ARGV = @argv;

my %seq;
while (<>) {
    next unless /\S/;
    next if /^\s*\#/;
    if (/^\s*\/\//) { printseq() }
    else {
	chomp;
	my ($name, $seq) = split;
	$seq =~ s/[\.\-]//g unless $gapped;
	$seq{$name} .= $seq;
    }
}
printseq();

sub printseq {
	if($sorted){
		foreach $key (sort keys %seq){
			print ">$key\n";
			for (my $i = 0; $i < length $seq{$key}; $i += $columns){
				print substr($seq{$key}, $i, $columns), "\n";
			}
		}
	} else{
    		while (my ($name, $seq) = each %seq) {
			print ">$name\n";
			for (my $i = 0; $i < length $seq; $i += $columns) {
	    			print substr ($seq, $i, $columns), "\n";
			}
   	 	}
	}
    %seq = ();
}
