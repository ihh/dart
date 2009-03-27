
=head1 NAME

PhyloGram.pm

=head1 SYNOPSIS

Perl module for working with xrate-format phylo-grammar files.

The module extends the DartSexpr.pm module, and so inherits all of the methods from that class.

For a description of the format itself, see the following URL:

http://biowiki.org/XrateFormat

=head1 EXAMPLES

=head2 1. Estimate a general irreversible substitution model from an alignment

 use Stockholm::Database;
 use PhyloGram::Dna;

 $d=PhyloGram::Dna->new;
 $c=$d->new_chain("A");
 $d->add_emission("A X");
 $d->add_transition("X","X");
 $d->add_end_transition("X");

 $s=Stockholm::Database->from_file("rfam/Training/Cobalamin.stock");
 $d->train($s);
 print $d->find_chain("A")->mutate("a","c")->rate->value,"\n";

=head2 2. Find conserved alignment regions using a 10-state phylo-HMM

 use Stockholm::Database;
 use PhyloGram::Dna;

 $d=PhyloGram::Dna->new;
 $d->new_chain("A");
 $d->add_emission("A X");
 $d->add_transition("X","X");
 $d->add_end_transition("X");
 $s=Stockholm::Database->from_file("rfam/Training/Lysine.stock");
 $d->xgram_args("-log ECFG_STATS");
 $d->train($s);
 $c=$d->find_chain("A");

 $g=PhyloGram::Dna->new;
 $g->add_end_transition("S")->prob(.001);
 foreach$k(1..10){
  $ck=$g->new_empty_chain("k$k");
  $ck->inc_initial($c);
  $ck->inc_rates($c,$k/10);
  $g->add_emission("k$k x$k")->add(["annotate",["row","Conservation"],["column","k$k"],["label",chr(48+10-$k)]]);
  $g->add_transition("S","x$k")->prob(.0999);
  $g->add_transition("x$k","S");
 }
 $g->annotate($s);
 print $s->to_string(-1),"\n";

=head2 3. Display type table of codon grammar

 use PhyloGram;
 $t=PhyloGram->from_file("$ENV{DARTDIR}/grammars/codon.eg")->symbol_type;
 foreach$s(sort{$$t{$a}cmp$$t{$b}||$a cmp$b}keys%$t){
  print"$s => $$t{$s}\n";
 }

=head2 4. Print out all the rate matrices in a grammar

 use PhyloGram;
 $g=PhyloGram->from_file("$ENV{DARTDIR}/grammars/codon.eg");
 foreach$chain($g->all_chains){
  print$chain->terminal->to_string,"\n";
  $mutateHash=$chain->mutate_hash;
  @states=map(join("",@{$_->state->value}),$chain->find_all("initial"));
  print"* @states\n";
  foreach$i(@states){
   print$i;
   foreach$j(@states){
    print" ",$i eq$j?"*":$chain->mutate($i,$j,$mutateHash)->rate->value;
   }
   print"\n";
  }
  print"-"x80,"\n";
 }

=head2 5. Print out all the parameter-value pairs in a grammar

 use PhyloGram;
 $g=PhyloGram->from_file("$ENV{DARTDIR}/grammars/hky85.eg");
 $h=$g->param_hash;
 while(($n,$s)=each%$h){
  print$n," ",$s->value,"\n";
 }

=head1 METHODS

=cut

package PhyloGram;

use DartSexpr;
use PhyloGram::Chain;
use Stockholm::Database;

@ISA = qw(DartSexpr);

use strict;
use vars '@ISA';

use Carp;

# xrate executables
my $xrateExec = "xrate";
my $dartlog = "dartlog.pl";

# special chars: post_emit_char & complement_char
sub post_emit_char { "*" }
sub obsolete_post_emit_char { "'" }  # for backward compatibility
sub complement_char { "~" }

my $COMPLEMENT_PATTERN = "\\" . complement_char;
my $POST_EMIT_PATTERN = "[\\"  . post_emit_char . "\\" . obsolete_post_emit_char . "]";

=head2 new_default_grammar

    my $gram = PhyloGram->new_default_grammar();

Creates an empty grammar.

=cut

# constructor
sub new_default_grammar {
    my ($class) = @_;

    # create S-expression
    my $self = DartSexpr->new;
    bless $self, $class;

    # add grammar
    $self->add (['grammar',
		 ['name' => 'PhyloGram-default-grammar']
		 ]);

    # return
    return $self;
}

=head2 tidy_string

    print $gram->tidy_string


This method formats the grammar neatly. Note that this is done by "training" on an empty database, using xrate.
Therefore some information may be dropped (i.e. anything xrate does not recognize).

=cut

# tidy_string method: formats the grammar by "training" on an empty database
sub tidy_string {
    my ($self) = @_;

    my ($grammarfile, $alignfile) = ($self->tmp_grammar_filename, $self->tmp_alignment_filename);
    $self->to_file ($grammarfile);
    Stockholm::Database->new->to_file ($alignfile);

    $self->run_xgram ("$alignfile --grammar $grammarfile --train $grammarfile --noannotate", 10);

    my $tidy_string = `cat $grammarfile`;

    unlink $grammarfile if -e $grammarfile;
    unlink $alignfile if -e $alignfile;

    return $tidy_string;
}

=head2 train

    $gram->train ($stockholmDatabase)
    $gram->train ($stockholmDatabase, $commandLineOptions)

This method "trains" the grammar on a Stockholm::Database object, by running xrate.
If any command-line options are specified, they will be passed to xrate.

=cut

# train methods: use 'xgram' to optimize parameters for a Stockholm database
sub train {
    my ($self, $stockdb, $cloptions) = @_;
    $cloptions = "" unless defined $cloptions;

    my ($grammarfile, $alignfile) = ($self->tmp_grammar_filename, $self->tmp_alignment_filename);
    $self->to_file ($grammarfile);
    $stockdb->to_file ($alignfile);

    $self->run_xgram ("$alignfile --grammar $grammarfile --train $grammarfile --noannotate $cloptions");

    my $newgram = DartSexpr->from_file ($grammarfile);
    @$self = @$newgram;

    unlink $grammarfile if -e $grammarfile;
    unlink $alignfile if -e $alignfile;

    return $self;
}

=head2 annotate

    $gram->annotate ($stockholmDatabase)
    $gram->annotate ($stockholmDatabase, $commandLineOptions)

This method "annotates" a Stockholm::Database object from the grammar, by running xrate.
If any command-line options are specified, they will be passed to xrate.

=cut

# annotate method: use 'xgram' to annotate a Stockholm database
sub annotate {
    my ($self, $stockdb, $cloptions) = @_;

    my ($grammarfile, $alignfile) = ($self->tmp_grammar_filename, $self->tmp_alignment_filename);
    $self->to_file ($grammarfile);
    $stockdb->to_file ($alignfile);

    my $newstock_text = $self->run_xgram ("$alignfile --grammar $grammarfile --annotate $cloptions");
    my $newstockdb = Stockholm::Database->from_string ($newstock_text);
    @$stockdb = @$newstockdb;

    unlink $grammarfile if -e $grammarfile;
    unlink $alignfile if -e $alignfile;

    return $stockdb;
}

=head2 estimate_missing_trees_by_nj

    $gram->estimate_missing_trees_by_nj ($stockholmDatabase)

This method uses neighbor-joining to estimate trees for a Stockholm::Database object, by running xrate.
The branch lengths of the trees are then optimized using the EM algorithm.

=cut

# estimate_missing_trees_by_nj method: use 'xgram' to add trees to a Stockholm database
# (topology and branch lengths of any missing trees will be estimated by neighbor-joining)
# warning: default mutation rates are of order .01, so branches will be long if you use default chain constructors
sub estimate_missing_trees_by_nj {
    my ($self, $stockdb) = @_;

    my ($grammarfile, $alignfile) = ($self->tmp_grammar_filename, $self->tmp_alignment_filename);
    $self->to_file ($grammarfile);
    $stockdb->to_file ($alignfile);

    my $newstock_text = $self->run_xgram ("$alignfile --tree $grammarfile --noannotate");
    my $newstockdb = Stockholm::Database->from_string ($newstock_text);
    @$stockdb = @$newstockdb;

    unlink $grammarfile if -e $grammarfile;
    unlink $alignfile if -e $alignfile;

    return $stockdb;
}

# estimate_branch_lengths_by_em method: use 'xgram' to re-estimate branch lengths of trees in a Stockholm database
# (topology of any missing trees will be estimated by neighbor-joining)
# warning: default mutation rates are of order .01, so branches will be long if you use default chain constructors
sub estimate_branch_lengths_by_em {
    my ($self, $stockdb) = @_;

    my ($grammarfile, $alignfile) = ($self->tmp_grammar_filename, $self->tmp_alignment_filename);
    $self->to_file ($grammarfile);
    $stockdb->to_file ($alignfile);

    my $newstock_text = $self->run_xgram ("$alignfile --branches --grammar $grammarfile --noannotate");
    my $newstockdb = Stockholm::Database->from_string ($newstock_text);
    @$stockdb = @$newstockdb;

    unlink $grammarfile if -e $grammarfile;
    unlink $alignfile if -e $alignfile;

    return $stockdb;
}

=head2 new_pgroup

    my $sexpr = $gram->new_pgroup (@param)

This method adds a new pgroup to the grammar.
The parameter names are given in the @param array.
The default value for each parameter is 1/@param.

The return value is the S-expression for the new pgroup.

=cut

# new_pgroup method
sub new_pgroup {
    my ($self, @param) = @_;
    return $self->grammar->params->add_child (map ([$_ => 1/@param], @param));
}

=head2 new_rate

    my $sexpr = $gram->new_rate ($rateParam, $value)

This method adds a new rate parameter to the grammar.

The return value is the S-expression for the new rate parameter.

=cut

# new_rate method
sub new_rate {
    my ($self, $param, $val) = @_;
    $val = 1 unless defined $val;
    return $self->grammar->params->add_child ($param => $val);
}

=head2 new_empty_chain

    my $chain = $gram->new_empty_chain (@pseudoterminals)

This method adds a new Markov chain with the specified pseudoterminals.

The return value is a PhyloGram::Chain object.

=cut

# new_empty_chain method
# all rates & probabilities zero
sub new_empty_chain {
    my ($self, @pseudoterm) = @_;
    croak "No pseudoterminal array" unless @pseudoterm;
    @pseudoterm = map (split, @pseudoterm);
    my $chain = DartSexpr->new ('chain',
				['update-policy' => 'irrev'],
				['terminal' => \@pseudoterm],
				);
    $chain = $self->bless_chain ($chain);
    $self->grammar->add ($chain);
    return $chain;
}

=head2 new_sparse_chain

    my $chain = $gram->new_sparse_chain (@pseudoterminals)

This method adds a new Markov chain with the specified pseudoterminals.
The chain is then "sparse-filled", i.e. for every single-residue mutation a default entry is added.

The return value is a PhyloGram::Chain object.

=cut

# new_sparse_chain method
# creates & sparse-fills a chain
sub new_sparse_chain {
    my ($self, @pseudoterm) = @_;
    my $chain = $self->new_empty_chain (@pseudoterm);
    $self->chain_init ($chain);
    $self->chain_sparse_fill ($chain);
    return $chain;
}


=head2 new_chain

    my $chain = $gram->new_chain (@pseudoterminals)

This method adds a new Markov chain with the specified pseudoterminals.
The chain is then "dense-filled", i.e. for every possible mutation (single- or multi-residue) a default entry is added.

The return value is a PhyloGram::Chain object.

=cut

# new_chain method
# creates & dense-fills a chain
sub new_chain {
    my ($self, @pseudoterm) = @_;
    my $chain = $self->new_empty_chain (@pseudoterm);
    $self->chain_init ($chain);
    $self->chain_fill ($chain);
    return $chain;
}

=head2 new_dense_chain

    my $chain = $gram->new_dense_chain (@pseudoterminals)

A synonym for the new_chain method.

=cut

# new_dense_chain = new_chain
sub new_dense_chain { new_chain (@_) }

=head2 param_hash

    my %param_hash = %{$gram->param_hash()}

This method returns a reference to a hash, mapping parameter names to S-expression nodes where they are defined.

The method is useful because otherwise, retrieving the value of a particular parameter may involve trawling through the entire grammar file,
which is extremely inefficient and slow.

=cut

# param_hash method
# returns a hashref of parameter names to S-expression nodes where they are defined
sub param_hash {
    my ($self) = @_;
    my @paramList =  # loop over children of grammar node
	$self->grammar->map_child (sub {
	    my ($grammar_child) = @_;
	    # check that this grammar-child is a list, whose first element is an atom
	    if (ref ($grammar_child) && !ref ($grammar_child->[0])) {
		# check if the atom is one of the parameter-declaration keywords
		my $tag = $grammar_child->[0];
		if ($tag =~ /^(params|const|pgroup|const-pgroup|rate|const-rate)$/) {
		    # loop over children of the parameter-declaration node
		    return $grammar_child->map_child (sub {
			my ($child) = @_;
			# return a list of parameter-value nodes
			# children of the parameter-declaration node come in two types:
			# pgroup nodes: ((A .1) (B .2) (C .3)  ... )
			#   rate nodes: (R 10)
			return
			    ref($child)
			    ? (@$child >= 2 && !ref($child->[0]) && !ref($child->[1])  # rate node?
			       ? ($child)  # rate (single)
			       : @$child)  # pgroup (list)
			    : ();
		    });
		} else {
		    return ();
		}
	    }
	});
    my %paramHash = map (($_->[0] => $_), @paramList);
    return \%paramHash;
}

=head2 find_chain

    my $chain = $gram->find_chain ($pseudoterminal)
    my @chains = $gram->find_chain ($pseudoterminal)

    my $codeRef = sub {
	my $chain = shift;
	...
	    # return true or false
	};

    my $chain = $gram->find_chain ($codeRef);
    my @chains = $gram->find_chain ($codeRef);

This method returns a chain (or list of chains) containing a given pseudoterminal, or matching a given predicate.

Each chain is returned as a PhyloGram::Chain object.

=cut

# find_chain method
# finds chain(s) with given pseudoterminal (or matching given predicate, if a CODE ref is supplied)
sub find_chain {
    my ($self, $pseudoterminal) = @_;
    my $codeRef = ref($pseudoterminal) eq 'CODE'
	? $pseudoterminal
	: sub {
	    my $chain = shift;
	    return 0 + grep ($_ eq $pseudoterminal, @{$chain->terminal->value});
	};
    my @chains = grep (&$codeRef($_), $self->all_chains);
    return wantarray() ? @chains : $chains[0];
}

=head2 all_chains

    my @chains = $gram->all_chains()

This method returns all the chains in the grammar.

Each chain is returned as a PhyloGram::Chain object.

=cut

# all_chains method
# returns an array of all chains
sub all_chains {
    my ($self) = @_;
    return map ($self->bless_chain($_), $self->grammar->find_all ('chain'));
}

# re-bless a chain... seems kind of hacky and unprincipled, so isolate it here
sub bless_chain {
    my ($self, $chain) = @_;
    bless $chain, 'PhyloGram::Chain' if ref($chain) ne 'PhyloGram::Chain';
    return $chain;
}

=head2 chain_states

    my @states = $gram->all_chains ($chain)

This method returns all the states for a given PhyloGram::Chain object.

=cut

# set of all states for a chain
sub chain_states {
    my ($self, $chain) = @_;

    my @tok = @{$self->alphabet->token->value};
    my $terms = @{$chain->terminal->value} + 0;

    my $ntok = @tok + 0;
    my $states = $ntok**$terms;

    my $int2state = sub {
	my ($int) = @_;
	return [map ($tok[int($int / ($ntok ** $_)) % $ntok], 0..$terms-1)];
    };
    my @states = map (&$int2state($_), 0..$states-1);

    return @states;
}

=head2 sparse_mutations

    my @destStates = $gram->sparse_mutations ($srcState, $position)

This method returns the list of all "sparse" (single-residue) mutations from a given state at a given position.

=cut

# return the list of all "sparse" mutations from a given state at given position
sub sparse_mutations {
    my ($self, $from, $pos) = @_;
    my @tok = @{$self->alphabet->token->value};
    my @to;
    foreach my $tok (@tok) {
	if ($tok ne $$from[$pos]) {
	    my $to = [@$from];
	    $$to[$pos] = $tok;
	    push @to, $to;
	}
    }
    return @to;
}

=head2 chain_init

    $gram->chain_init ($chain)

This method sets all initial probabilities for a PhyloGram::Chain object to (1/total_number_of_states).

=cut

# chain_init method
# sets all initial probabilities to 1.0 / (number of states)
sub chain_init {
    my ($self, $chain) = @_;
    $chain = $self->bless_chain ($chain);  # just in case

    my @states = $self->chain_states ($chain);
    foreach my $state (@states) {
	$chain->initial ($state, 1/@states);
    }
}

=head2 chain_sparse_fill

    $gram->chain_sparse_fill ($chain, $rate)

This method sparse-fills a chain. See new_sparse_chain() method for details.

=cut

# chain_sparse_fill method
# initialises rates for all single-residue mutations in a chain to a const value
sub chain_sparse_fill {
    my ($self, $chain, $rate) = @_;
    $chain = $self->bless_chain ($chain);  # just in case

    my @states = $self->chain_states ($chain);
    $rate = .01/@states unless defined $rate;   # default rate is .01/states (start small => faster convergence)

    my $mutate_hash = $chain->mutate_hash;
    foreach my $from (@states) {
	for (my $pos = 0; $pos < @$from; ++$pos) {
	    my @to = $self->sparse_mutations ($from, $pos);
	    foreach my $to (@to) {
		$chain->mutate ($from, $to, $rate, $mutate_hash);
	    }
	}
    }
}

=head2 chain_fill

    $gram->chain_fill ($chain, $rate)

This method dense-fills a chain. See new_dense_chain() method for details.

=cut

# chain_fill method
# initialises rates for all mutations in a chain to a const value
# NB for single-terminal chains, chain_sparse_fill and chain_fill are identical
sub chain_fill {
    my ($self, $chain, $rate) = @_;

    my @states = $self->chain_states ($chain);
    $rate = .01/@states unless defined $rate;  # default rate is .01/states (start small => faster convergence)

    my $mutate_hash = $chain->mutate_hash;
    for (my $i = 0; $i < @states; ++$i) {
	for (my $j = 0; $j < @states; ++$j) {
	    if ($j != $i) {
		$chain->mutate ($states[$i], $states[$j], $rate, $mutate_hash);
	    }
	}
    }
}

=head2 chain_dense_fill

    $gram->chain_dense_fill ($chain, $rate)

This method dense-fills a chain. It is a synonym for the chain_fill method.

=cut

# chain_dense_fill = chain_fill
sub chain_dense_fill { chain_fill (@_) }

=head2 add_emission

    my $transformationRule = $gram->add_emission ($rhs)

This method adds a new transformation rule to the grammar, corresponding to an emission.

The $rhs argument should correspond to the right-hand side of the rule,
and consist of a whitespace-delimited string of pseudoterminals plus one nonterminal.

=cut

# add_emission
sub add_emission {
    my ($self, $rhs) = @_;
    $rhs =~ s/([\(\)])/ $1 /g;
    my @rhs = split (/\s+/, $rhs);
    grep (s/\s//g, @rhs);

    my ($lhs, $chain, @lcontext, @emit, @rcontext);
    my $rhsArrayRef = \@emit;
    foreach my $rhsSymbol (@rhs) {

	if ($rhsSymbol eq '(') {
	    @lcontext = @emit;
	    @emit = ();

	} elsif ($rhsSymbol eq ')') {
	    $rhsArrayRef = \@rcontext;

	} else {

	    # strip off complement & post-emit characters
	    # TODO: remove hardcoded chars here by using $self->post_emit_char and $self->complement_char
	    my $strippedRhsSymbol = $rhsSymbol;
	    $strippedRhsSymbol =~ s/^$COMPLEMENT_PATTERN//;  # strip off leading complement character, if it's there
	    $strippedRhsSymbol =~ s/[$POST_EMIT_PATTERN]$//;  # strip off trailing post-emit character, if it's there

	    # find chain
	    my $rhsChain = $self->find_chain ($strippedRhsSymbol);
	    if (defined $rhsChain) {
		croak "Chain mismatch on RHS of emit rule" if defined($chain) && $rhsChain != $chain;
		$chain = $rhsChain;

	    } else {
		croak "Multiple nonterminals ($lhs,$strippedRhsSymbol) on RHS of emit rule" if defined($lhs) && $strippedRhsSymbol ne $lhs;
		$lhs = $strippedRhsSymbol;
		$rhsSymbol = $lhs . $self->post_emit_char;
	    }

	    push @$rhsArrayRef, $rhsSymbol;
	}
    }

    croak "No LHS symbol" unless defined $lhs;

    # add emit rule
    my $transform = DartSexpr->new ('transform',
				    ['from' => [@lcontext, $lhs, @rcontext]],
				    ['to' => [@lcontext, @emit, @rcontext]],
				    );

    $self->grammar->add ($transform);
    return $transform;
}

=head2 emission

    my $emission = $gram->emission ($nonterminal)

Returns the transformation rule corresponding to the emission from a given nonterminal, or undef.

=cut

# method to test for emit state
sub emission {
    my ($self, $state) = @_;

    my $post_emit = $state . $self->post_emit_char;
    my @emission = $self->grammar->grep_child (sub {
	my ($child) = @_;
	return $child->has_tag && $child->tag eq 'transform'
	    && grep ($_ eq $state, @{$child->from->value})
	    && grep ($_ eq $post_emit, @{$child->to->value});
    });

    die map ($_->to_string . "\n", @emission), "Multiple emissions for state $state" if @emission > 1;
    return @emission ? $emission[0] : undef;
}

=head2 add_transition

    my $transition = $gram->add_transition ($from, $to, $prob)

Adds a transition to the grammar.

=cut

# add_transition
sub add_transition {
    my ($self, $from, $to, $prob) = @_;
    $prob = "1" unless defined $prob;
    $from .= $self->post_emit_char if defined $self->emission ($from);

    # add transition
    my $transform = DartSexpr->new ('transform',
				    ['from' => [$from]],
				    ['to' => defined($to) ? [$to] : []],
				    ['prob' => $prob],
				    );

    $self->grammar->add ($transform);
    return $transform;
}

=head2 add_end_transition

    my $end_transition = $gram->add_transition ($from)

Adds an end transition to the grammar.

=cut

# add_end_transition
sub add_end_transition {
    my ($self, $from) = @_;
    return $self->add_transition ($from, undef);
}

=head2 add_bifurcation

    my $bifurc = $gram->add_bifurcation ($from, $toLeft, $toRight)

Adds a bifurcation to the grammar.

=cut

# add_bifurcation
sub add_bifurcation {
    my ($self, $from, $to_l, $to_r) = @_;
    die "Can't bifurcate from emit state $from" if defined $self->emission ($from);

    # add bifurcation
    my $transform = DartSexpr->new ('transform',
				    ['from' => [$from]],
				    ['to' => [$to_l, $to_r]],
				    );

    $self->grammar->add ($transform);
    return $transform;
}

=head2 transition

    my $transition = $gram->transition ($from, $to)

Locates a transition between two given states, adding it to the grammar if it was not already there.

=cut

# transition accessor
# calls add_transition if not found
sub transition {
    my ($self, $from, $to) = @_;
    $from .= $self->post_emit_char if $self->emission ($from);

    my @trans = $self->grammar->grep_child (sub {
	my ($child) = @_;
	return $child->has_tag && $child->tag eq 'transform'
	    && @{$child->from->value} == 1 && $child->from->value->[0] eq $from
	    && (defined($to)
	    ? (@{$child->to->value} == 1 && $child->to->value->[0] eq $to)
	    : (@{$child->to->value} == 0));
    });

    croak "Multiple transform clauses from $from to $to" if @trans > 1;
    unless (@trans) {
	push @trans, $self->add_transition ($from, $to);
    }

    return $trans[0];
}

=head2 end_transition

    my $endTransition = $gram->end_transition ($from)

Locates a transition from a given state to the end state, adding it to the grammar if it was not already there.

=cut

# end_transition accessor
# calls add_transition if not found
sub end_transition {
    my ($self, $from) = @_;
    return $self->transition ($from, undef);
}

=head2 symbol_type

    my %type = $gram->symbol_type()

Returns the symbol type table of the grammar: a map from identifiers to the following strings:

 alphabet
 token
 extend
 terminal
 emit
 null

=cut

# symbol type table
sub symbol_type {
    my ($self) = @_;

    my $comp_re = "^" . $COMPLEMENT_PATTERN . "(.+)\$";
    my $post_re = "^(.+)" . $POST_EMIT_PATTERN . "\$";

    my %alphabet = map (($_ => 'alphabet'),
			map ($_->name->value,
			     $self->find_all ('alphabet')));

    my %token = map (($_ => 'token'),
		     map (@{$_->token->value},
			  $self->find_all ('alphabet')));

    my %extend = map (($_ => 'extend'),
		      map ($_->to->value,
			   map ($_->find_all ('extend'),
				$self->find_all ('alphabet'))));

    my %terminal = map (($_ => 'terminal'),
			map (@{$_->terminal->value},
			     map ($_->find_all ('chain'),
				  $self->find_all ('grammar'))));

    my %emit = map (($_ => 'emit'),
		    grep (!exists ($terminal{$_}) && s/$post_re/$1/,
			  grep (s/$comp_re/$1/ || 1,
			  map (@{$_->to->value},
			       map ($_->find_all ('transform'),
				    $self->find_all ('grammar'))))));

    my %null = map (($_ => 'null'),
		    grep (!exists ($terminal{$_}) && !exists ($emit{$_}),
			  grep (s/$comp_re/$1/ || 1,
				grep (s/$post_re/$1/ || 1,
				      map ((@{$_->from->value}, @{$_->to->value}),
					   map ($_->find_all ('transform'),
						$self->find_all ('grammar')))))));

    return {%alphabet, %token, %extend, %terminal, %emit, %null};
}

=head2 new_symbol

    my $newSymbolName = $gram->new_symbol()
    my $newSymbolName = $gram->new_symbol($preferredName)
    my $newSymbolName = $gram->new_symbol($preferredName,$type)

Method to generate new unique symbol.

=cut

# method to generate new unique symbol
sub new_symbol {
    my ($self, $preferred, $symbol_type, $new_type) = @_;
    $preferred = "S" unless defined $preferred;
    $symbol_type = $self->symbol_type unless defined $symbol_type;

    my $name = $preferred;
    for (my $n = 2; exists $symbol_type->{$name}; ++$n) {
	$name = $preferred . $n;
    }

    $symbol_type->{$name} = $new_type;
    return $name;
}

=head2 tmp_filename

    my $tmpFilename = $gram->tmp_filename ($prefix)

Method to generate new unique temporary filename.

=head2 tmp_grammar_filename

    my $tmpFilename = $gram->tmp_grammar_filename()

Does what you might expect.

=head2 tmp_alignment_filename

    my $tmpFilename = $gram->tmp_alignment_filename()

Does what you might expect.

=head2 tmp_log_filename

    my $tmpFilename = $gram->tmp_log_filename()

Does what you might expect.

=cut

# temporary filename method
sub tmp_filename {
    my ($prefix) = @_;
    my $suffix = $$;
    my $filename = $prefix . $suffix;
    for (my $inc = 1; -e $filename; ++$inc) {
	$filename = "$prefix$suffix.$inc";
    }
    return $filename;
}

# temporary grammar, alignment & log filenames
sub tmp_grammar_filename { return tmp_filename ("phylo-grammar.") }
sub tmp_alignment_filename { return tmp_filename ("phylo-aligndb.") }
sub tmp_log_filename { return tmp_filename ("phylo-log.") } # return "xgram.log" }

=head2 run_xrate

    my $output = $gram->run_xrate ($arguments, $logLevel)

Runs xrate with the specified arguments and logging level.

=cut

# method to run xgram
sub run_xgram {
    my ($self, $args, $loglevel) = @_;
    $loglevel = 6 unless defined $loglevel;   # default log level; lower for more info

    my $logfile = $self->tmp_log_filename;
    my $existing_logfile = -e $logfile;

    my @xgram_args = map ($_->values, $self->find_all ('xgram-args'));
    my $command = "$xrateExec $args -log $loglevel -logfile $logfile @xgram_args";
    warn "...running '$command'\n";

    my ($output, $child);
    if ($child = fork) {

	$SIG{INT} = $SIG{KILL} = sub { kill 'KILL', $child };

	$output = `$command`;
	kill (9, $child);

	delete $SIG{INT};
	delete $SIG{KILL};

    } else {

	# dartlog.pl call commented out because it seems to leave zombies -- IH, 1/25/2008
	# TODO: fix this

#	while (!-e $logfile) { }  # sleep(1)

#	local *DARTLOG;
#	local $_;
#	open DARTLOG, "$dartlog -tail $logfile |";
#	$SIG{INT} = $SIG{KILL} = sub { close DARTLOG };

#	while (<DARTLOG>) { print STDERR $_ }

#	close DARTLOG;
#	$SIG{INT} = $SIG{KILL} = undef;

	exit;
    }

    unlink $logfile if -e $logfile && !$existing_logfile;

    return $output;
}

# synonym
sub run_xrate {
    return run_xgram (@_);
}

# cool builder methods

=head2 hmm_emit

    my ($nonterminal, $chain_sexpr, $emit_sexpr) = $gram->hmm_emit()
    my ($nonterminal, $chain_sexpr, $emit_sexpr) = $gram->hmm_emit ($preferredNonterminal)
    my ($nonterminal, $chain_sexpr, $emit_sexpr) = $gram->hmm_emit ($preferredNonterminal, $pseudoterminal)
    my ($nonterminal, $chain_sexpr, $emit_sexpr) = $gram->hmm_emit ($preferredNonterminal, $pseudoterminal, $update_policy)

Adds a left-emit state to a grammar and automatically creates an associated chain.

=cut

# method to add a left-emit state and automatically create an associated chain
# assumes a single grammar.
sub hmm_emit {
    my ($self, $nonterm, $pseudoterm, $update_policy) = @_;
    my $symbol_type = $self->symbol_type;

    $nonterm = $self->new_symbol (defined($nonterm) ? $nonterm : "hmm-emit", $symbol_type, 'emit');
    $pseudoterm = $self->new_symbol (defined ($pseudoterm) ? $pseudoterm : "\$" . $nonterm, $symbol_type, 'terminal');
    $update_policy = 'rind' unless defined $update_policy;

    my $chain_sexpr = $self->new_chain ($pseudoterm);
    $chain_sexpr->update_policy ($update_policy);
    my $emit_sexpr = $self->add_emission ("$pseudoterm $nonterm");

    # return state, chain S-expr, transformation rule S-expr
    return ($nonterm, $chain_sexpr, $emit_sexpr);
}

=head2 ghmm_emit

    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit()
    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit ($size)
    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit ($size, $firstNonterminal)
    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit ($size, $firstNonterminal, $pseudoterminal)
    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit ($size, $firstNonterminal, $pseudoterminal, $update_policy)
    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit ($size, $firstNonterminal, $pseudoterminal, $update_policy, $emitNonterminal)
    my ($firstNonterminal, $pseudoterminal, $lastNonterminal, @emitNonterminals) = $gram->hhmm_emit ($size, $firstNonterminal, $pseudoterminal, $update_policy, $emitNonterminal, $endNonterminal)

Adds a whole range of emit states gatewayed by null states, allowing for length distributions (effectively a Generalised HMM).

=cut

# method to add a whole range of emit states gatewayed by null states, allowing for length distributions (effectively a Generalised HMM)
sub ghmm_emit {
    my ($self, $size, $start_nonterm, $pseudoterm, $update_policy, $emit_nonterm, $end_nonterm) = @_;
    my $symbol_type = $self->symbol_type;

    $size = 1000 unless defined $size;   # default is 1kb
    $start_nonterm = $self->new_symbol (defined($start_nonterm) ? $start_nonterm : "ghmm", $symbol_type, 'null');
    $pseudoterm = $self->new_symbol (defined($pseudoterm) ? $pseudoterm : "\$" . $start_nonterm, $symbol_type, 'terminal');
    $update_policy = 'rev' unless defined $update_policy;
    $emit_nonterm = defined($emit_nonterm) ? $emit_nonterm : "$start_nonterm-emit";
    $end_nonterm = $self->new_symbol (defined($end_nonterm) ? $end_nonterm : "$start_nonterm-end", $symbol_type, 'null');

    # create chain
    my $chain_sexpr = $self->new_chain ($pseudoterm);
    $chain_sexpr->update_policy ($update_policy);

    # create emit states
    my $prev_nonterm = $end_nonterm;
    my @emit_nonterm;
    for (my $i = 0; $i < $size; ++$i) {
	my $current_nonterm = $self->new_symbol ($emit_nonterm,
						 $symbol_type, 'emit');
	$self->add_emission ("$pseudoterm $current_nonterm");
	$self->add_transition ($start_nonterm, $current_nonterm)->prob (1 / $size);
	$self->add_transition ($current_nonterm, $prev_nonterm);

	push @emit_nonterm, $current_nonterm;
	$prev_nonterm = $current_nonterm;
    }

    # return state & chain ids
    return ($start_nonterm, $pseudoterm, $end_nonterm, @emit_nonterm);
}

# end

1;
