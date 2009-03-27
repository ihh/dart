#!/usr/bin/perl

use strict;
use PhyloGram::Protein;
use Stockholm::Database;

# Number of rate classes
my $N = 10;

# Label of annotation row
my $annotRow = "ANNOT";

# Read in the nullprot grammar
my $nullprotFile = "../nullprot.eg";  # hardcode for now
my $nullprot = PhyloGram->from_file ($nullprotFile);

# Extract the rates from the nullprot grammar
my ($nullprot_chain) = $nullprot->all_chains;
my $nullprot_mutate_hash = $nullprot_chain->mutate_hash;  # this is sufficient for later random-access to rates
my @states = map (join ("",@{$_->state->value}), $nullprot_chain->find_all("initial"));  # get all nullprot states

# Create protein grammar object
my $hmm = PhyloGram::Protein->new;
$hmm->grammar->name ("$N-nonterminal-protein-grammar");
$hmm->grammar->parametric (1);

# Add the nullprot initial probs & rates to $hmm as constant parameters
my @constPi;
for my $state (@states) {
    push @constPi, [piParam($state),$nullprot_chain->initial($state)->prob->value];
}
my @constParam = (\@constPi);
for my $src (@states) {
    for my $dest (@states) {
	if ($dest ne $src) {
	    my $mutate = $nullprot_chain->mutate($src,$dest,$nullprot_mutate_hash);
	    push @constParam, [[rateParam($src,$dest), defined($mutate) ? $mutate->rate->value : 0]];
	}
    }
}
$hmm->grammar->add (["const",@constParam]);

# Add the N scaling rates, and the mixture weights, to $hmm as variable parameters
my @param = map ([[scaleParam($_),$_/$N]],1..$N);
for my $k (1..$N) {
    push @param, [[weightParam($k),1/$N]];
}
$hmm->grammar->add (["params",@param]);

# Add transitions from Start to Null, and Null to End
$hmm->add_transition ("Start", "Null", 1);
$hmm->add_transition ("Null", "", 1);


# Create the N emit nonterminals; populate the parametric rate matrices with the appropriate scaled rate functions
for my $k (1..$N) {
    warn "Creating nonterminal $k\n";
    # declare the k'th emit nonterm & the k'th chain
    my ($nonterm, $chain, $rule) = $hmm->hmm_emit (nontermName($k), termName($k), "parametric");
    # add transitions to/from Null
    # (this must happen AFTER the emit state is declared...)
    $hmm->add_transition ("Null", nontermName($k), weightParam($k));
    $hmm->add_transition (nontermName($k), "Null", 1);
    # Add annotation clause
    $rule->add (["annotate",["row",$annotRow],["column",termName($k)],["label",columnLabel($k)]]);
    # set initial probabilities
    for my $state (@states) {
	$chain->initial ($state, piParam($state));
    }
    # set mutation rates
    my $mutate_hash = $chain->mutate_hash;
    for my $src (@states) {
	for my $dest (@states) {
	    if ($dest ne $src) {
		$chain->mutate ($src, $dest, scaleParam($k) . " * " . rateParam($src,$dest), $mutate_hash);
	    }
	}
    }
}

# Output the grammar
warn "Formatting grammar file\n";
print $hmm->tidy_string;

# parameter, nonterminal, terminal names; labels
sub rateParam { my ($src, $dest) = @_; return "r" . uc($src) . uc($dest) }
sub piParam { my ($state) = @_; return "pi" . uc($state) }
sub weightParam { my ($k) = @_; return "w$k" }
sub scaleParam { my ($k) = @_; return "s$k" }
sub nontermName { my ($k) = @_; return "Nonterm$k" }
sub termName { my ($k) = @_; return "Term$k" }
sub columnLabel { my ($k) = @_; return chr($k+64) }
