#!/usr/bin/perl -w

package ChainDat;

use strict;
use vars '@ISA';

use Carp;

use PhyloGram;

########## Handy variables ##########
my %nuc = (
	   A => 1,
	   C => 2,
	   G => 3,
	   T => 4,
	  );
my %dinuc = (
	     AA => 1,
	     AC => 2,
	     AG => 3,
	     AT => 4,
	     CA => 5,
	     CC => 6,
	     CG => 7,
	     CT => 8,
	     GA => 9,
	     GC => 10,
	     GG => 11,
	     GT => 12,
	     TA => 13,
	     TC => 14,
	     TG => 15,
	     TT => 16,
	    );
my %codon = (  # don't handle stop codons
	     AAA => 1,
	     AAC => 2,
	     AAG => 3,
	     AAT => 4,
	     ACA => 5,
	     ACC => 6,
	     ACG => 7,
	     ACT => 8,
	     AGA => 9,
	     AGC => 10,
	     AGG => 11,
	     AGT => 12,
	     ATA => 13,
	     ATC => 14,
	     ATG => 15,
	     ATT => 16,
	     CAA => 17,
	     CAC => 18,
	     CAG => 19,
	     CAT => 20,
	     CCA => 21,
	     CCC => 22,
	     CCG => 23,
	     CCT => 24,
	     CGA => 25,
	     CGC => 26,
	     CGG => 27,
	     CGT => 28,
	     CTA => 29,
	     CTC => 30,
	     CTG => 31,
	     CTT => 32,
	     GAA => 33,
	     GAC => 34,
	     GAG => 35,
	     GAT => 36,
	     GCA => 37,
	     GCC => 38,
	     GCG => 39,
	     GCT => 40,
	     GGA => 41,
	     GGC => 42,
	     GGG => 43,
	     GGT => 44,
	     GTA => 45,
	     GTC => 46,
	     GTG => 47,
	     GTT => 48,
	     TAC => 49,
	     TAT => 50,
	     TCA => 51,
	     TCC => 52,
	     TCG => 53,
	     TCT => 54,
	     TGC => 55,
	     TGG => 56,
	     TGT => 57,
	     TTA => 58,
	     TTC => 59,
	     TTG => 60,
	     TTT => 61,
	    );
my %aa = (
	  A => 1,
	  C => 2,
	  D => 3,
	  E => 4,
	  F => 5,
	  G => 6,
	  H => 7,
	  I => 8,
	  K => 9,
	  L => 10,
	  M => 11,
	  N => 12,
	  P => 13,
	  Q => 14,
	  R => 15,
	  S => 16,
	  T => 17,
	  V => 18,
	  W => 19,
	  Y => 20,
	  );


# constructor
sub new {
    my ($class, $params, $chain) = @_;
    my $self = {
		'params' => {}, # parameter hashref
		'chain' => [],
		'N' => "",
	       };
    bless $self, $class;
    $self->params($params);
    $self->chain($chain);
    return $self->_initialize;
}

sub _initialize {
  my ($self) = @_;

  my @states = map(join("", @{$_->state->value}), $self->chain->find_all("initial"));
  $self->N(scalar(@states)); # compute alphabet size

  return $self;
}

# catch methods
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD; # $AUTOLOAD contains the fully qualified name of the original subroutine
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq 'DESTROY';

  # check for params accessor, e.g. $self->params_Ka or $self->params_('Ka')
  if ($sub =~ /^params_(\S*)$/i) {
    my $param = lc($1);
    $param = shift @args unless (length $param); # this catches the argument 'Ka' in the second example usage given above
    if (!defined $self->{"params"}->{$param}) {
      $self->{"params"}->{$param} = ""; # if no such parameter exists, create one with the value ""
    }
    return $self->{"params"}->{$param};
  }

  # check for ordinary accessors
  # This has the effect of automatically implementing getter and setter methods.
  # If there's an instance variable $name, then the getter method is $self->name
  # and the setter method is $self->name('newName')
  if (exists $self->{$sub}) {
    if (@args > 1) { croak "Usage: $sub() or $sub(newValue)"; }
    return
      @args                      # if @args > 0
      ? $self->{$sub} = $args[0] # { $self->{$sub} = $args[0]; return $args[0]; }
      : $self->{$sub};           # else { return $self->{$sub}; }
  }

  croak "Unsupported method: $sub\n";
}


sub initToDat {
  my ($self, $dat) = @_;
  my $N = $self->N;

  my @states = map(join("", @{$_->state->value}), $self->chain->find_all("initial"));
  my $stopCodonFlag = 0;
  foreach my $i (@states) {
    my $newi = uc($i); # for simplicity, use uppercase DNA alphabet
    $newi =~ s/U/T/ig;
    if    ($N == 4) { print $dat "$nuc{$newi} "; }
    elsif ($N == 16) { print $dat "$dinuc{$newi} "; } # don't handle stop codons
    elsif (($N == 61) || ($N == 64)) { if (($newi eq "TAA") || ($newi eq "TAG") || ($newi eq "TGA")) { $stopCodonFlag = 1; next; }
		       print $dat "$codon{$newi} "; }
    elsif ($N == 20) { print $dat "$aa{$newi} "; }
    if ($self->chain->initial($i)) {
      my @values = $self->chain->initial($i)->prob->values;
      my $value = $self->parseParams(\@values);
      print $dat "$value\n";
    } else { print $dat "0\n"; } # init. prob. not specified
  }
  if ($stopCodonFlag) { print "\n\nWARNING: I'm ignoring the stop codons in your grammar!\n\n\n"; }
}


sub ratesToDat {
  my ($self, $dat) = @_;
  my $N = $self->N;

  # for fast indexing, it's necessary to make a lookup table of the rate matrix. The PhyloGram::Chain method "mutate_hash" does this
  my $mutateHash = $self->chain->mutate_hash;
  # get a list of states by extracting the "state" field of all child nodes named "initial" ("find_all" and "value" are DartSexpr methods)
  my @states = map(join("", @{$_->state->value}), $self->chain->find_all("initial"));
  foreach my $i (@states) {
    my $newi = uc($i);
    $newi =~ s/U/T/ig;
    foreach my $j (@states) {
      unless ($i eq $j) {
	my $newj = uc($j); # for simplicity, use uppercase DNA alphabet
	$newj =~ s/U/T/ig;
	if    ($N == 4)                  { print $dat "$nuc{$newi} $nuc{$newj} "; }
	elsif ($N == 16)                 { print $dat "$dinuc{$newi} $dinuc{$newj} "; } # don't handle stop codons
	elsif (($N == 61) || ($N == 64)) { if (($newi eq "TAA") || ($newi eq "TAG") || ($newi eq "TGA")) { next; }
					   if (($newj eq "TAA") || ($newj eq "TAG") || ($newj eq "TGA")) { next; }
					   print $dat "$codon{$newi} $codon{$newj} "; }
	elsif ($N == 20)                 { print $dat "$aa{$newi} $aa{$newj} "; }
	if ($self->chain->mutate($i, $j, $mutateHash)) { # if the mutation rate is specified in the grammar
	  my @values = $self->chain->mutate($i, $j, $mutateHash)->rate->values;
	  my $value = $self->parseParams(\@values);
	  print $dat "$value\n";
	} else { print $dat "0\n"; } # rate not specified
      }
    }
  }

}


sub parseParams {
  my ($self, $values) = @_;
  my $value;
  
  if (@$values > 1) {		# if mathematical expression...
    for (my $k = 0; $k < @$values; $k++) {
      # is it a reference? if so, call recursively
      if (ref ($values->[$k])) {
	$values->[$k] = $self->parseParams($values->[$k]);
      }
      # is it a defined parameter? if so, substitute
      if ($self->params_($values->[$k])) {
	$values->[$k] = $self->params_($values->[$k]);
      }
      # unless decimal or scientific notation number or mathematical operator, it's illegal
      unless ($values->[$k] =~ /^\d*\.?\d*$/
	      || $values->[$k] =~ /^\d*\.?\d*e-\d+$/
	      || $values->[$k] eq "*"
	      || $values->[$k] eq "/"
	      || $values->[$k] eq "+"
	      || $values->[$k] eq "-") {
	croak "\n\nThe parameter '$values->[$k]' in '",@$values,"' isn't defined in your grammar!\n\n\n";
      }
    }
    $value = eval(join("",@$values));
  }
  else {			# single parameter or value
    $value = $values->[0];

    # check for special case of entire expression enclosed in parentheses
    # is it a reference? if so, call recursively
    if (ref ($value)) {
      $value = $self->parseParams($value);
    }

    # is it a defined parameter?  if so, substitute
    if ($self->params_($value)) {
      $value = $self->params_($value);
    }
    # unless decimal or scientific notation number, it's illegal
    unless ($value =~ /^\d*\.?\d*$/
	    || $value =~ /^\d*\.?\d*e-\d+$/) {
      croak "\n\nThe parameter '$value' isn't defined in your grammar!\n\n\n";
    }
  }

  return $value;
}

1
