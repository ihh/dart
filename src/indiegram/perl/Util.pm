#!/usr/bin/perl -w

use strict;

package Util;


################
# Math methods #
################

# Returns min (a, b).
sub min {
  my ($class, $a, $b) = @_;
  if ($a < $b) { return $a; }
  return $b;
}

# Log base 2.
sub logTwo {
  my ($class, $n) = @_;
  return (log ($n) / log (2.));
}


##################
# String methods #
##################


# Strip out all non-alphanumeric characters from a string.
sub toAlphanumeric {
  my ($class, $string) = @_;
  $string =~ s/\W//g;
  return $string;
}

# Convert a string into an array of characters; return reference to said array
sub toArray {
  my $str = shift;

  my @arr = split(//,$str);
  return \@arr;
}

# safer version of substr (doesn't complain if startpoint outside of string)
# (written by Ian)
sub safe_substr {
    my ($string, $start, $len) = @_;
    return "" if !defined($string) || $start > length $string;
    return substr ($string, $start, $len);
}

#################
# Array methods #
#################


# Given two array refs, computes the intersection of the 
# arrays and returns the corresponding array ref.
sub isect {
  my ($class, $a, $b) = @_;

  my %union;
  my %isect;
  foreach my $n (@$a, @$b) {
    $union{$n}++ && $isect{$n}++;
  }

  my @isect = keys %isect;
  return \@isect;
}


# Given two array refs, computes the union of the 
# arrays and returns the corresponding array ref.
sub union {
  my ($class, $a, $b) = @_;

  my %union;
  foreach my $n (@$a, @$b) {
    $union{$n}++;
  }

  my @union = keys %union;
  return \@union;
}


# Given an array ref and a member, return true if the array
# contains the member and false if it doesn't.
sub aContains {
  my ($class, $aRef, $n) = @_;

  foreach my $m (@$aRef) {
    if ($m eq $n) { return 1; }
  }
  return 0;
}

# Given an array ref and a member of the array, find and remove each 
# entry which is equal to the member.
sub aRemove {
  my ($class, $aRef, $n) = @_;

  # find all instances of $n
  my @indices;
  for (my $i = 0; $i < @$aRef; $i += 1) {
    if ($n eq $aRef->[$i]) { push (@indices, $i); }
  }
  
  # remove them from the array
  foreach my $index (@indices) {
    splice (@$aRef, $index, 1);
  }
}


sub deepCopy {
  my ($class, $this) = @_;
  if (not ref $this) {
    $this;
  } elsif (ref $this eq "ARRAY") {
    [map $class->deepCopy ($_), @$this];
  } elsif (ref $this eq "HASH") {
    +{map { $_ => $class->deepCopy ($this->{$_}) } keys %$this};
  } else { die "'$_' isn't 'ARRAY' or 'HASH'; I'm bailing.\n" }
}

1
