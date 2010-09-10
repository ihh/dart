#!/usr/bin/env perl -w

use strict;
use Stockholm;

my $usage = "\nUsage: $0 <ClustalW file>
           Converts a ClustalW file (from file or STDIN) to Stockholm format.
           Only converts sequence information.\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-/) {
	if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
	else { die $usage; }
    } else {
	push @argv, $arg;
    }
}
push @argv, "-" unless @argv;
my $file = shift @argv or die $usage;

open ALN, "<$file" or die "Couldn't open '$file'\n";
my $stk = Stockholm->new();

my $gapChars = '-._';
while (<ALN>) {
  my @a = split;
  next unless @a == 2;

  my ($seq, $data) = @a;
  next if ($data =~ /[^a-zA-Z$gapChars]/); # skip primary sequence conservation lines, etc.

  if (defined $stk->seqdata->{$seq}) {
    $stk->seqdata->{$seq} .= $data;
  } else {
    push @{$stk->seqname}, $seq;
    $stk->seqdata->{$seq} = $data;
  }

}
die unless $stk->is_flush();

print $stk->to_string();
