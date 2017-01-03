#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use strict;
use Stockholm;

my $usage = "\nUsage: $0 <column format (FoldAlign) file>
           Converts a column format (FoldAlign) file (from file or STDIN) to Stockholm format.
           Converts sequence and consensus structure information (forces <> notation); ignores all else.\n";

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

open COL, "<$file" or die "Couldn't open '$file'\n";
my $stk = Stockholm->new();
# initialize SS_cons line
$stk->gc->{"SS_cons"} = "";

# Sequence data is in the format:
# ; ALIGN               X06054        GGGCCCGUCG UCUAGCCUGG UUAG-GACGC UGCCCUGACG

my $gapChars = '-._';
my %seqs; # $seqs{"seq"} defined means that this seq is present in the alignment
while (<COL>) {
  next unless /^;\sALIGN\s/;
  chomp;

  my @a = split;

  # skip empty lines
  if (@a <= 2) { next; }

  # use the Begin line to get sequence names
  if (/Begin/i) {
    my $seq = $a[2];
    push (@{$stk->seqname}, $seq);
    $stk->seqdata->{$seq} = ""; # initialize sequence data
    $seqs{$seq} = 1;
    next;
  }
  # skip the End line
  elsif (/End/i) {
    next;
  }
  # are we looking at the consensus structure line?
  elsif ($a[2] eq "Structure") {
    for (my $i = 3; $i < @a; $i++) {
      # convert from () to <>
      my $ss = $a[$i];
      $ss =~ tr/\(\)/<>/;
      $stk->gc->{"SS_cons"} .= $ss;
    }
  }
  # are we looking at sequence data?
  elsif (defined $seqs{$a[2]}) {
    my $seq = $a[2];
    for (my $i = 3; $i < @a; $i++) {
      $stk->seqdata->{$seq} .= $a[$i];
    }
  }
  else {
    warn "Skipping:\n$_";
    next;
  }
}

die unless $stk->is_flush();
print $stk->to_string();
