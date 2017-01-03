#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use PhyloGram;
use Carp;

use strict;

my $usage = "\nUsage: $0 [options] <Xrate grammar file>

           [-h, --help] display this message
           [--const] force parameters to be constant (default is not)
           [--symmetric] force the chain to be strand-symmetric (default is not)

This script takes a (optionally partially parameterized) grammar file as input, 
extracts all probabilities and then stores them as parameters in
'pgroup' and 'rate' blocks.  If requested, it will also make
the chains strand-symmetric (so that e.g. 'A' and 'T' have identical initial probabilities
and the rates 'A -> T' and 'T -> A' are identical).

Writes to STDOUT.

NB: The \"running 'xrate ...'\" warning can be ignored.  This is a hack for formatting the grammar.
\n";

my $infile;
# set defaults
my $const = 0;
my $symm = 0;

while (@ARGV) {
  my $arg = shift @ARGV;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif ($arg eq "--const") { $const = 1; }
    elsif (($arg eq "--symmetric")) { $symm = 1; }
    else { croak $usage; }
  }
  $infile = $arg;
}
unless (defined $infile) { print $usage; exit; }

# load grammar
my $gram = PhyloGram->from_file ($infile);

# get (original) parameters
my %oldparams;
my $h = $gram->param_hash;
while (my ($n,$s) = each %$h) { $oldparams{$n} = $s->value;}

# if it's an RNA alphabet, use U instead of T
my $rna = 0;
if ($gram->alphabet->name->value =~ /rna/i) { $rna = 1; }


foreach my $chain ($gram->all_chains()) {

# read in init. prob. dist. and rate matrix for each chain of the grammar
my %init;
my %rate;

# used for symmetrization
my %init_new;
my %rate_new;

  $chain->tag_value("update-policy", "parametric");
  # get states of the chain
  my @states = map (join ("", @{$_->state->value}), $chain->find_all("initial"));

  # for fast indexing, it's necessary to make a lookup table of the rate matrix. The PhyloGram::Chain method "mutate_hash" does this
  my $mutateHash = $chain->mutate_hash;

  # set prefix for init and rate probs to the pseudoterminal symbols of the chain
  my $prefix = join ("_",@{$chain->terminal->[1]});


  # get the initial dist
  for my $i (@states) {
    my $probname = "$prefix"."_$i";
    my @thisprob = $chain->initial($i)->prob->values;
    my $parsedprob = parseParams(\%oldparams, \@thisprob);
    $init{$probname} = $parsedprob;

    # get the rates
    for my $j (@states) {
      next if $i eq $j;
      if (!defined $chain->mutate ($i,$j,$mutateHash)) { warn "missing $i -> $j!\n"; next; }
      my @thisrate = $chain->mutate($i,$j,$mutateHash)->rate->values;
      my $parsedrate = parseParams(\%oldparams, \@thisrate);
      my $ratename = "$prefix"."_$i"."2"."$j";
      $rate{$ratename} = $parsedrate;
    }
  }


  # Now tag the grammar and symmetrize as we go if requested.
  # If not symmetrizing, then copy %init into %init_new.
  if (!$symm) {
    while (my ($n,$p) = each %init) { $init_new{$n} = $p; }
    while (my ($n,$p) = each %rate) { $rate_new{$n} = $p; }
  }
  
  # loop through the states, tagging the grammar as we go
  for my $i (@states) {
    my $probname = "$prefix"."_$i";
    # if not symmetrizing, then just tag the grammar
    if (!$symm) {
      $chain->initial($i)->tag_value("prob", $probname);
    }
    
    # else do lots of funny stuff to enforce symmetrization
    else {
      my $probname_wc = $prefix . "_" . complement ($i);
      my $newprobname = $probname . "_" . complement ($i);
      my $newprobname_wc = $probname_wc . "_" . $i;

      # if we've already calculated the init. prob. (as its WC equivalent), then just tag the grammar
      if (defined $init_new{$newprobname_wc}) {
#	$chain->initial($i)->tag_value("prob", $newprobname_wc);
	$chain->initial($i)->tag_value("prob", "($newprobname_wc / 2)");
      } 
      # else store the init. prob. and tag the grammar
      else {
	$init_new{$newprobname} = ($init{$probname} + $init{$probname_wc});
	$chain->initial($i)->tag_value("prob", "($newprobname / 2)");
      }
    }

    for my $j (@states) {
      next if $i eq $j;
      if (!defined $chain->mutate ($i,$j,$mutateHash)) { warn "missing $i -> $j!\n"; next; }
      my $ratename = "$prefix"."_$i"."2"."$j";
      # if not symmetrizing, then just tag the grammar
      if (!$symm) {
	$chain->mutate($i,$j,$mutateHash)->tag_value("rate", $ratename);
      }

      # else enforce symmetrization
      else {
	my $ratename_wc = $prefix . "_" . complement ($i."2".$j);
	my $newratename = $ratename . "_" . complement ($i."2".$j);
	my $newratename_wc = $ratename_wc . "_" . $i."2".$j;

	# if we've already calculated the rate (as its WC equivalent), then just tag the grammar
	if (defined $rate_new{$newratename_wc}) {
	  $chain->mutate($i,$j,$mutateHash)->tag_value("rate", $newratename_wc);
	}
	# else store the rate and tag the grammar
	else {
	  $rate_new{$newratename} = ($rate{$ratename} + $rate{$ratename_wc}) / 2;
	  $chain->mutate($i,$j,$mutateHash)->tag_value("rate", $newratename);
	}
      }
    }

  }

# modify the grammar to add in the new parameters
	my @finalinits;
	my $initSexpr;
	my $i=0;
	foreach my $key (sort keys %init_new) {
		if($symm){
#			my $thisProb = $init_new{$key} / 2; 
#			$finalinits[$i] = "(".$key." ".$thisProb.")";
			$finalinits[$i] = "(".$key." ".$init_new{$key}.")";
		}
		else{
  			$finalinits[$i] = "(".$key." ".$init_new{$key}.")";
		}
  		$i++;
	}		
	foreach my $sexper (@finalinits) {
  		$initSexpr .= $sexper;
	}

	if($const) {
 		my $constpgroup = "const-pgroup";
  		$gram->grammar->find_or_add($constpgroup);
  		$gram->grammar->$constpgroup->add_child($initSexpr);
	} else {
  		$gram->grammar->find_or_add("pgroup");
  		$gram->grammar->pgroup->add_child($initSexpr);
	}

	if ($const){
  		my $constrate = "const-rate";
  		$gram->grammar->find_or_add($constrate);
  		foreach my $key (sort keys %rate_new){
    			$gram->grammar->$constrate->add_child($key." ".$rate_new{$key});
  		}
	}else{
  		$gram->grammar->find_or_add("rate");
  		foreach my $key (sort keys %rate_new){
    			$gram->grammar->rate->add_child($key." ".$rate_new{$key});
 	 	}
	}
} # end loop over chains





# print the result to STDOUT
print $gram->tidy_string();



##########################subs#################################################################

# Return the Watson-Crick complement of the passed string.
sub complement {
  my ($s) = @_;
  if ($rna) { $s =~ tr/acgtu/ugcaa/; }
  else { $s =~ tr/acgtu/tgcaa/; }
  $s = reverse $s;
  return $s;
}


#modified from ChainDat.pm - LB
sub parseParams {
  my ($params, $values) = @_;
  my $value;
  
  if (@$values > 1) {		# if mathematical expression...
    for (my $k = 0; $k < @$values; $k++) {
      # is it a reference? if so, call recursively
      if (ref ($values->[$k])) {
		$values->[$k] = parseParams($params, $values->[$k]);
      }
      # is it a defined parameter? if so, substitute
      if ($params->{$values->[$k]}) {
		$values->[$k] = $params->{$values->[$k]};
      }
      # unless decimal or scientific notation number or mathematical operator, it's illegal
      unless ($values->[$k] =~ /^\d*\.?\d*$/g
	      || $values->[$k] =~ /^\d*\.?\d*e-\d+$/g
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
      $value = parseParams($params,$value);
    }

    # is it a defined parameter?  if so, substitute
    if ($params->{$value}) {
      $value = $params->{$value};
    }
    # unless decimal or scientific notation number, it's illegal
    unless ($value =~ /^\d*\.?\d*$/g
	    || $value =~ /^\d*\.?\d*e-\d+$/g) {
      croak "\n\nThe parameter '$value' isn't defined in your grammar!\n\n\n";
    }
  }

  return $value;
}
