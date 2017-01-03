#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use strict;

use Stockholm;
use Stockholm::Database;

my $usage = "\nUsage: $0 <Stockholm alignment>
         [-h, --help] display this message
         [--pairwise] report average pairwise percent ID instead of default
         [--countgaps] use gap information (see below)
         [--foreach] print percent ID for each alignment (default is to average over alignments)
Calculates per-column % id as: (# of times most common non-gap character in column appears) / (# non-gap characters in column).
Only counts columns with > 1 non-gap character unless 'countgaps' flag is set.
If 'countgaps' flag is set, becomes (# of times most common non-gap character in column appears) / (# total characters in column).
Beware: the average percent ID calculation takes a long time!\n";

my $pairwise = 0;
my $countgaps = 0;
my $foreach = 0;

my @argv;
while (@ARGV) {
  my $arg = shift;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif (($arg eq "--pairwise")) { $pairwise = 1; }
    elsif (($arg eq "--countgaps")) { $countgaps = 1; }
    elsif (($arg eq "--foreach")) { $foreach = 1; }
    else { die $usage; }
  } else {
    push @argv, $arg;
  }
}

push @argv, "-" unless @argv;
my $file = shift @argv or die $usage;
warn "[reading alignments]\n";
my $db = Stockholm::Database->from_file ($file);

my $sumID = 0; # sum of percent ID's defined over single alignments
my $n = 0; # number of percent ID calculations

warn "[processing ", @$db+0, " alignments]\n";
for my $stock (@$db) {

  if (!$stock->is_flush()) { die "Not flush; I'm quitting\n"; }

  my $cols = $stock->columns;
  warn "[processing $cols columns]\n";

  if ($pairwise) {
    for (my $s1 = 0; $s1 < @{$stock->seqname}; $s1++) {
      for (my $s2 = $s1+1; $s2 < @{$stock->seqname}; $s2++) {
	my ($seqname1, $seqname2) = ($stock->seqname->[$s1], $stock->seqname->[$s2]);

	warn "[processing $seqname1 $seqname2]\n";
	my @seqdata = ($stock->seqdata->{$seqname1}, $stock->seqdata->{$seqname2});
	my $id = percentID (\@seqdata);
	$sumID += $id;
	++$n;
      }
    }
  }

  else {
    my @seqdata;
    while (my ($name,$seq) = each %{$stock->seqdata}) { push @seqdata, $seq; }
    my $id = percentID (\@seqdata);
    $sumID += $id;
    ++$n;
  }

  if ($foreach) {
    print "Identity: ", 100 * $sumID / $n, " percent\n";
    $sumID = 0;
    $n = 0;
  }

}

if (!$foreach) {
  print "Identity: ", 100 * $sumID / $n, " percent\n";
}


# Calculate percentID for the passed set of sequences.
# Calculates per-column % id as: (# of times most common non-gap character in column appears) / (# non-gap characters in column)
sub percentID {
  my ($seqs) = @_;
  my $cols = length $seqs->[0];
  my $rows = scalar @$seqs;

  my $sum = 0;
  my $n = 0;

  for (my $c = 0; $c < $cols; ++$c) {
    my %f;
    my $count = 0; # ungapped
    my $countall = 0; # gaps as well
    for (my $i = 0; $i < $rows; ++$i) {
      my $char = lc substr ($seqs->[$i], $c, 1);
      if ($char ne '-' && $char ne '.') {
	++$f{$char};
	++$count;
      }
      ++$countall;
    }
    # unless $countgaps, only count columns with > 1 non-gap character
    if ((!$countgaps) && ($count > 1)) {
      my @sym = sort { $f{$b} <=> $f{$a} } keys %f;

      my $id;
      # if no character appears more than once, then % id = 0
      if ($f{$sym[0]} == 1) {
	$id = 0;
      }
      # else calculate % id as (# most common character) / (# non-gap characters in column)
      else {
	$id = $f{$sym[0]} / $count;
      }
      $sum += $id;
      ++$n;
    }
    # if $countsgaps, count columns with >= 1 non-gap character
    elsif (($countgaps) && ($count >= 1)) {
      my @sym = sort { $f{$b} <=> $f{$a} } keys %f;

      my $id;
      # if no character appears more than once, then % id = 0
      if ($f{$sym[0]} == 1) {
	$id = 0;
      }
      # else calculate % id as (# most common character) / (# total characters in column)
      else {
	$id = $f{$sym[0]} / $countall;
      }
      $sum += $id;
      ++$n;
    }

  }

  if ($n == 0) {
    die "Insufficient data";
  }
  
  return $sum / $n;
}
