#!/usr/bin/perl -w

package TreeTransducer::BranchTT;

use Carp;
use strict;
use TreeTransducer;

our @ISA = qw /TreeTransducer/;
use vars '@ISA';

sub new {
  my ($class, $filename) = @_;

  my $self = TreeTransducer->from_file ($filename);

  bless $self, ref ($class) || $class;

  $self->_checkValid();

  return $self;
}


# Various checks that the branch transducer is valid.
sub _checkValid {
  my ($self) = @_;

  my %validTypes = (
		    's' => 1,
		    'i' => 1,
		    'm' => 1,
		    'bi' => 1, # insert bifurcation
		    'bm' => 1, # match bifurcationg
		    'w' => 1,
		    'e' => 1
		   );
  
  # Check that every state has a valid state type.
  foreach my $s (@{$self->states()}) {
    my $t = $self->stateType ($s);
    if (!$t || !exists $validTypes{$t}) {
	croak ("Empty or invalid state typing found for (branch) state '$s'.\n");
    }
  }

}


1
