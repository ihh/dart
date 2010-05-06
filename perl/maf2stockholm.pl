#!/usr/bin/env perl -w

use Stockholm;

# params
my $cols = 50;

# help message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname [-c <column width>] <MAF file>\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-c") {
	defined ($cols = shift) or die $usage;
    } elsif ($arg =~ /^-/) {
	die $usage, "\nUnrecognized argument: $arg\n";
    } else {
	push @argv, $arg;
    }
}
unless (@argv) {
    @argv = ('-');
    warn "[waiting for alignments on standard input]\n";
}
@ARGV = @argv;

my $stock;
while (<>) {
    next if /^\s*\#/;  # comment
    if (/^a /) {
	print $stock->to_string($cols) if defined $stock;
	$stock = Stockholm->new;
    } elsif (/^s /) {
	my ($s_char, $name, $start, $len, $strand, $end, $seq) = split;
	$stock->add_row ($name, $seq);
    }
}
print $stock->to_string($cols) if defined $stock;
