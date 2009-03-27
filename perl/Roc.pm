#!/usr/bin/perl -w

package Roc;

use strict;
use vars '@ISA';
use Carp;

# constructor
sub new {
    my $class = shift(@_);
    my $self = {
      cat => undef,
      tp => 0,
      tn => 0,
      fp => 0,
      fn => 0,
      nullCols => 0,
      @_ # Override previous attributes
    };
    bless $self, $class;
    return $self;
}

sub evalColPred
{
  my $self = shift(@_);
  my $actual = shift(@_);
  my $pred = shift(@_);
  if ($actual eq ".")
  {
    # Ignore these columns
    $self->{nullCols}++;
    return;
  }
  else
  {
    if ($actual eq $self->{cat})
    {
      if ($actual eq $pred) {$self->{tp}++}
      else {$self->{fn}++}
    }
    elsif ($pred eq $self->{cat}) {$self->{fp}++}
    else {$self->{tn}++}
  }
}

sub getNullCols
{
  my $self = shift(@_);
  return $self->{nullCols};
}

sub toString
{
  my $self = shift(@_);
  my $sn = $self->{tp}/($self->{tp} + $self->{fn}) * 100;
  my $sp = $self->{tn}/($self->{fp} + $self->{tn}) * 100;
  print "cat=$self->{cat}, tp=$self->{tp}, fn=$self->{fn}, tn=$self->{tn}, fp=$self->{fp}, "
    . "Sn=$sn%, Sp=$sp%\n";
}

1
