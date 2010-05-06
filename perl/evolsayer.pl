#!/usr/bin/env perl -w

# evolsayer: simulator for TKF Structure Tree
# Ian Holmes, 2008-12-08

# Known issues:
# 1. Slows down hugely for long sequences (despite optimization). Inevitable really -- long sequences have lots of mutations.
# 2. Random number seeding seems fragile, for some reason. (Changing the command-line args, or modifying the code, seems to give different random numbers.)
# 3. Initial TKFST generation is recursive, so Perl gives "Deep recursion..." warnings for big trees. Could fix this with an iterative implementation.

# standard library imports
use Carp;

# DART imports
use Newick;
use PhyloGram;
use Stockholm;

# special TKFST states
my ($loopState, $stemState, $deadState, $unbornState) = qw(L S D I);
my $specialStateRegexp = "[$loopState$stemState$deadState$unbornState]";

# secondary structure character constants
my ($lchar, $rchar, $sschar, $gapchar) = qw(< > . -);

# get program name & directory
my $progname = $0;
$progname =~ s!.*/!!;

my $progdir = $0;
$progdir =~ s!/[^/]+$!!;

# get DARTDIR and grammar dir
my $dartdir = $ENV{'DARTDIR'};
if (!defined $dartdir) {
    if ($progdir =~ /\//) {
	$dartdir = $progdir;
	$dartdir =~ s!/[^/]+$!!;
    } else {
	$dartdir = "$progdir/..";
    }
}

# filenames & suffices
my $historySuffix = ".history";

# initialize params (as of 12/8/2008, these are the evoldoer defaults)
my $params = Params->new;
$params->loop_length (9);
$params->mu_loop (.03);
$params->stem_length (2.3);  # evoldoer default is actually 7/3 = 2.333333333
$params->mu_stem (.01);
$params->stem_prob (.01);
$params->subfile ("$dartdir/grammars/pfold-mix80.eg");

# initialize other options
my $verbose = 0;
my $history = 0;
my $rndseed = time;
my $s_root = 0;

# initialise help string
my $usage = "$progname: simulator for the TKF Structure Tree.\n";
$usage .= "\n";
$usage .= "Usage:\n";
$usage .= "       $progname [options]  <Newick tree file>\n";
$usage .= "\n";
$usage .= "Options:\n";
$usage .= "\t              [-h]  this help message\n";
$usage .= "\t              [-v]  print verbose status messages during simulated evolution\n";
$usage .= "\t     [-seed <int>]  seed for random number generator (default is clock time)\n";
$usage .= "\t [-sub <filename>]  substitution model filename (default is ".$params->subfile.")\n";
$usage .= "\t[-looplen <float>]  mean length of loops (default ".$params->loop_length.")\n";
$usage .= "\t[-loopdel <float>]  deletion rate in loops (default ".$params->mu_loop.")\n";
$usage .= "\t[-stemlen <float>]  mean length of stems (default ".$params->stem_length.")\n";
$usage .= "\t[-stemdel <float>]  deletion rate in stems (default ".$params->mu_stem.")\n";
$usage .= "\t  [-pstem <float>]  probability of finding a stem inside a loop (default ".$params->stem_prob.")\n";
$usage .= "\t          [-sroot]  force at least one stem, by rooting TKFST in state S rather than state L\n";
$usage .= "\t        [-history]  record .history files during simulated evolution\n";
$usage .= "\n";
$usage .= "For description of TKF Structure Tree model, see:\n";
$usage .= "\n";
$usage .= "  Holmes I, 2004. A probabilistic model for the evolution of RNA structure.\n";
$usage .= "  BMC Bioinformatics 5:166.   http://www.biomedcentral.com/1471-2105/5/166\n";
$usage .= "\n";

# parse cmd-line opts
my @argv;
my @oldARGV = @ARGV;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-v") {
	$verbose = 1;
    } elsif ($arg eq "-seed") {
	defined ($rndseed = shift) or die $usage;
    } elsif ($arg eq "-sub") {
	defined ($params->subfile (shift)) or die $usage;
    } elsif ($arg eq "-looplen") {
	defined ($params->loop_length (shift)) or die $usage;
    } elsif ($arg eq "-loopdel") {
	defined ($params->mu_loop (shift)) or die $usage;
    } elsif ($arg eq "-stemlen") {
	defined ($params->stem_length (shift)) or die $usage;
    } elsif ($arg eq "-stemdel") {
	defined ($params->mu_stem (shift)) or die $usage;
    } elsif ($arg eq "-pstem") {
	defined ($params->stem_prob (shift)) or die $usage;
    } elsif ($arg eq "-sroot") {
	$s_root = 1;
    } elsif ($arg eq "-history") {
	$history = 1;
    } else {
	push @argv, $arg;
    }
}

unless (@argv) {
    @argv = ('-');
    warn "[waiting for tree on standard input]\n";
}

die $usage if @argv != 1;
my ($treefile) = @argv;

# seed random number generator & log to stderr
Node::srandprob ($rndseed);

# load tree
my $tree = Newick->from_file ($treefile);

# load substitution models
$params->read_rates();

# create TKFST
warn "Sampling root sequence\n" if $verbose;
my $tkfst = $s_root ? Node::sample_equilibrium_stem ($params) : Node::sample_equilibrium_loop ($params);

# record number of ancestral stems
my ($n_ancestral_stems, $n_nonempty_ancestral_stems) = (0, 0);
$tkfst->traverse (sub { if ($_[0]->state eq 'S') { ++$n_ancestral_stems; ++$n_nonempty_ancestral_stems if $_[0]->child->[0]->state ne 'L' } });

# record lengths of ancestral stems & loops
my ($anc_loop_len, $anc_stem_len) = get_loop_and_stem_lengths ($tkfst);

# evolve TKFST
my @tag_summary = $tkfst->evolve_tree ($params, $tree);

# create Stockholm alignment
my $stock = $tkfst->make_stockholm ($tree);

# add params to alignment
$stock->add_gf ("ARGS", join (" ", @oldARGV));
$stock->add_gf ("SEED", $rndseed);
$stock->add_gf ("PARAM", join (" ", map ("$_=".eval("\$params->$_"), qw(stem_prob subfile loop_length stem_length))));
$stock->add_gf ("PARAM", join (" ", map ("$_=".eval("\$params->$_"), qw(lambda_loop mu_loop lambda_stem mu_stem))));

# add event log to alignment
for my $tag_summary (@tag_summary) {
    $stock->add_gf (@$tag_summary);
}

# add other helpful info to alignment
$stock->add_gf ("ANCSTEMS", "$n_ancestral_stems (total) $n_nonempty_ancestral_stems (nonempty)");
$stock->add_gf ("ANCSTEMLEN", join (" ", @$anc_stem_len));
$stock->add_gf ("ANCLOOPLEN", join (" ", @$anc_loop_len));

# print
print $stock->to_string;

# exit
exit;

# helper method to count loop & stem lengths
sub get_loop_and_stem_lengths {
    my ($tkfst) = @_;
    my (%seqlen, %seqroot);
    my %typelen = map (($_ => []), qw(L S));
    $tkfst->traverse
	(sub {
	    my ($node) = @_;
	    if ($node->state =~ /^([LS])$/) {
		$seqroot{$node} = $node;
		$seqlen{$node} = 0;
	    } else {
		$seqroot{$node} = $seqroot{$node->parent};
		++$seqlen{$seqroot{$node}} unless $node->state =~ /D/;
	    }
	 },
	 sub { $seqroot{$_[0]} = $seqroot{$_[0]->parent} if defined $_[0]->parent },   # set S's root to parent L once stem has been processed
	 sub {
	     my ($node) = @_;
	     if ($node->state =~ /^([LS])$/) {
		 push @{$typelen{$1}}, $seqlen{$node};
	     }
	 });
    return map ($typelen{$_}, qw(L S));
}

# Params package
package Params;

sub new {
    my ($class) = @_;
    $class = ref($class) if ref($class);
    my $self = { 'subrate' => {},
		 'init_loop' => {},
		 'init_stem' => {} };
    bless $self, $class;
    return $self;
}

# lambdas
sub lambda_loop {
    my ($self) = @_;
    return $self->kappa_loop * $self->mu_loop;
}

sub lambda_stem {
    my ($self) = @_;
    return $self->kappa_stem * $self->mu_stem;
}

# kappas
sub kappa_loop {
    my ($self) = @_;
    return 1 - 1 / ($self->loop_length + 1);
}

sub kappa_stem {
    my ($self) = @_;
    return 1 - 1 / ($self->stem_length + 1);
}

# AUTOLOAD
sub AUTOLOAD {
    my ($self, @args) = @_;
    my $sub = our $AUTOLOAD;
    $sub =~ s/.*:://;  # strip off module path

    # check for DESTROY
    return if $sub eq "DESTROY";

    # accessor
    if (@args == 1) {
	$self->{$sub} = shift @args;
    } elsif (@args > 1) {
	$self->{$sub} = \@args;
    }
    return $self->{$sub};
}

# read_rates
# parses an xrate-format grammar file ($self->subfile), looking for two chains:
# a single-base rate matrix (1 pseudoterminal) and a base-pair rate matrix (2 pseudoterminals).
# (the PFOLD grammar, of course, has exactly this.)
# copies the rates into $self->subrate,
# and the initial probabilities into $self->init_stem and $self->init_loop.
# An extra character (S: probability $self->stem_prob) is added to the initial distribution for loops.
sub read_rates {
    my ($self) = @_;
    my $gram = PhyloGram->from_file ($self->subfile);
    for my $chain ($gram->all_chains) {
	my @term = @{$chain->terminal->values};
	my ($initref, $is_loop) = @term > 1 ? ($self->init_stem, 0) : ($self->init_loop, 1);
	my $mutateHash = $chain->mutate_hash;
	my @states;
	for my $init ($chain->find_all ("initial")) {
	    my $i = join ("", @{$init->state->value});
	    my $prob = ($is_loop ? 1 - $self->stem_prob : 1) * $init->prob->value;
	    $initref->{$i} = $prob;
	    if ($i =~ /$specialStateRegexp/) {
		die "Bad state: $i includes $specialStateRegexp";
	    }
	    push @states, $i;
	}
	if ($is_loop) {
	    $initref->{$stemState} = $self->stem_prob;
	}
	for my $i (@states) {
	    for my $j (@states) {
		if ($i ne $j) {
		    my $rate = 0;
		    my $mutate_sexpr = $chain->mutate ($i, $j, $mutateHash);
		    if (defined $mutate_sexpr) {
			$rate = $mutate_sexpr->rate->value;
		    }
		    $self->subrate->{$i}->{$j} = $rate;
		}
	    }
	}
    }
}

1;

# Node package
# represents a TKF structure tree as a doubly-linked tree.
package Node;

# constructor
sub new {
    my ($class, %hash) = @_;
    $class = ref($class) if ref($class);
    my $self = { 'parent' => undef,  # doubly-linked tree parent
		 'child' => [],  # doubly-linked tree children
		 'state' => undef,  # current TKFST state associated with this node
		 'cols' => 0,  # records whether initial state was zero-columns (L/S), one-column (a/c/g/u) or two-column (au/cg/gc/ua/etc)
		 'species_state' => {} };  # ancestral states
    while (my ($member, $value) = each %hash) {
	$self->{$member} = $value;
    }
    bless $self, $class;
    return $self;
}

# AUTOLOAD
sub AUTOLOAD {
    my ($self, @args) = @_;
    my $sub = our $AUTOLOAD;
    $sub =~ s/.*:://;  # strip off module path

    # check for DESTROY
    return if $sub eq "DESTROY";

    # accessor
    if (@args == 1) {
	$self->{$sub} = shift @args;
    } elsif (@args > 1) {
	$self->{$sub} = \@args;
    }
    return $self->{$sub};
}

# depth (used for debugging output)
sub depth {
    my ($node) = @_;
    my $depth = 0;
    while (defined $node->parent) {
	$node = $node->parent;
	++$depth if $node->state eq $loopState || $node->state eq $stemState;
    }
    return $depth;
}

# sample a tree structure from the equilibrium distribution
# currently split across several recursively-linked functions. TODO: make this iterative
sub sample_equilibrium_loop {
    my ($params, $parent) = @_;
#    warn "sample_equilibrium_loop: depth = ", $parent->depth if defined $parent;
    my $root = Node->new ("state" => $loopState, "parent" => $parent);
    push @{$parent->child}, $root if defined $parent;
    my $node = $root;
    while (randprob() < $params->kappa_loop) {
	$node = sample_loop_node ($params, $node);
    }
    return $root;
}

sub sample_equilibrium_stem {
    my ($params, $root) = @_;
#    warn "sample_equilibrium_stem: depth = ", $root->depth if defined $root;
    $root = Node->new ("state" => $stemState) unless defined $root;
    my $node = $root;
    while (randprob() < $params->kappa_stem) {
	$node = sample_stem_node ($params, $node);
    }
    sample_equilibrium_loop ($params, $node);
    return $root;
}

sub sample_loop_node {
    my ($params, $parent) = @_;
    my $node = Node->new (sample_state ($params->init_loop),
			  "parent" => $parent);
    push @{$parent->child}, $node if defined $parent;
#    warn "sample_loop_node: state = ", $node->state;
    if ($node->state eq $stemState) {
	sample_equilibrium_stem ($params, $node);
    }
    return $node;
}

sub sample_stem_node {
    my ($params, $parent) = @_;
    my $node = Node->new (sample_state ($params->init_stem),
			  "parent" => $parent);
#    warn "sample_stem_node: state = ", $node->state;
    push @{$parent->child}, $node if defined $parent;
    return $node;
}

sub sample_state {
    my ($initref) = @_;
    my $p = randprob();
    my $pcum = 0;
    my @state = keys %$initref;
    my $state;
    for my $s (@state) {
	my $prob = $initref->{$s};
	if (($pcum += $prob) >= $p) {
	    $state = $s;
	    last;
	}
    }
    if (!defined $state) {
	warn "sample_state fell out of loop (p=$p, total=$pcum); returning last state ($state[$#state])...";
	$state = $state[$#state];
    }
    my $cols = $state =~ /$specialStateRegexp/ ? 0 : length ($state);
    return ("state" => $state,
	    "cols" => $cols);
}

# random number stuff
sub srandprob {
    my ($seed) = @_;
    warn "Random number seed: $seed\n";
    srand ($seed);
}

sub randprob {
    return rand(1);
}

# tree traversal with callbacks for preorder (before visiting any children), mid-order (after each child) & postorder (after all children)
sub traverse {
    my ($node, $presubref, $midsubref, $postsubref) = @_;

    # here is the recursive code we are going to handle iteratively

    # &$presubref ($node) if defined $presubref;
    # my @child = @{$node->child};
    # for (my $n_child = 0; $n_child < @child; ++$n_child) {
    #	 $child[$n_child]->traverse ($presubref, $midsubref, $postsubref);
    #	 &$midsubref ($node, $n_child) if defined $midsubref;
    # }
    # &$postsubref ($node) if defined $postsubref;

    # here is the iterative version:
    my (@n_child_stack, @node_stack);
    &$presubref ($node) if defined $presubref;
    my $n_child = 0;
    while (defined $node) {
	if ($n_child < @{$node->child}) {
	    push @node_stack, $node;
	    push @n_child_stack, $n_child;
	    $node = $node->child->[$n_child];
	    $n_child = 0;
	    &$presubref ($node) if defined $presubref;
	} else {
	    &$postsubref ($node) if defined $postsubref;
	    $node = pop @node_stack;  # returns undef if node_stack is empty
	    $n_child = pop @n_child_stack;
	    if (defined $node) {
		&$midsubref ($node, $n_child) if defined $midsubref;
		++$n_child;
	    }
	}
    }
}

# methods to save and restore state, for simulation on a phylogenetic tree
# these store/retrieve the state of the entire TKFST in the species_state hashref,
# which is keyed by phylogenetic node ID.
sub save_state {
    my ($self, $species) = @_;
    $self->traverse
	(sub {
	    my ($node) = @_;
	    $node->species_state->{$species} = $node->state;
	 });
}

sub restore_state {
    my ($self, $species) = @_;
    $self->traverse
	(sub {
	    my ($node) = @_;
	    if (exists $node->species_state->{$species}) {
		$node->state ($node->species_state->{$species});
	    } else {
		$node->state ($unbornState);
	    }
	 });
}

# method to test if state is live
sub is_live {
    my ($state) = @_;
    return $state ne $deadState && $state ne $unbornState;
}

# method to copy entire TKFST structure, removing all non-live states
sub copy_live {
    my ($self) = @_;
    my @tkfst;  # new root(s)
    my $parent;  # current parent
    my @parent;  # parent stack
    $self->traverse (
	sub {
	    # pre
	    my ($node) = @_;
	    if (is_live ($node->state)) {
		my $copy = Node->new ("parent" => $parent,
				      "state" => $node->state,
				      "cols" => $node->cols,
				      "species_state" => { %{$node->species_state} });
		push @{$parent->child}, $copy if defined $parent;  # double-link the new tree
		push @tkfst, $copy unless defined $parent;  # if no parent, then this is a root node
		push @parent, $parent;  # put parent on stack
		$parent = $copy;
	    }
	},
	undef,
	sub {
	    # post
	    my ($node) = @_;
	    if (is_live ($node->state)) {
		$parent = pop @parent;  # get parent from stack
	    }
	});
    return @tkfst==1 ? $tkfst[0] : @tkfst;
}

# method to get alignment symbol for state, or undef
sub align_symbol {
    my ($state, $pos) = @_;
    return is_live ($state) ? substr ($state, $pos, 1) : undef;
}

# make_stockholm method
sub make_stockholm {
    my ($self, $tree) = @_;

    # get node IDs
    my @node_id = get_node_id ($tree);

    # create column lookups
    my (@col2seq, @col2ss, %species_col2seq);
    my (@lpos, @rpos, @rstack);
    my $col = 0;
    $self->traverse
	(sub {
	    # before visiting any children
	    my ($node) = @_;
	    my $state = $node->state;
#	    warn "state=$state depth=", $node->depth, " lpos=(@lpos) rpos=(@rpos)";
	    if ($state eq $stemState) {
		push @rstack, [@rpos];
		@rpos = ();
	    } elsif ($state ne $loopState) {
		if ($node->cols == 1) {
		    # unpaired column
		    push @lpos, $col;
		    $col2ss[$col] = $sschar;  #  "."
		    $col2seq[$col] = align_symbol ($state, 0);
		    while (my ($species, $state) = each %{$node->species_state}) {
			$species_col2seq{$species}->[$col] = align_symbol ($state, 0);
		    }
		    ++$col;

		} elsif ($node->cols == 2) {
		    # left column of pair
		    push @lpos, $col;
		    $col2ss[$col] = $lchar;  # "<"
		    $col2seq[$col] = align_symbol ($state, 0);
		    while (my ($species, $state) = each %{$node->species_state}) {
			$species_col2seq{$species}->[$col] = align_symbol ($state, 0);
		    }
		    ++$col;

		    # right column of pair
		    unshift @rpos, $col;
		    $col2ss[$col] = $rchar;  # ">"
		    $col2seq[$col] = align_symbol ($state, 1);
		    while (my ($species, $state) = each %{$node->species_state}) {
			$species_col2seq{$species}->[$col] = align_symbol ($state, 1);
		    }
		    ++$col;
		}
	    }
	 },
	 sub {
	    # after visiting each child
	    my ($node, $n_child) = @_;
	    my $state = $node->state;
	    if ($state eq $stemState && $n_child == 0) {
		push @lpos, @rpos;
		@rpos = @{pop @rstack};
#		warn "(n_child=$n_child) state=$state depth=", $node->depth, " lpos=(@lpos) rpos=(@rpos)";
	    }
	 });

    # get final order of columns
    my @cols = (@lpos, @rpos);

    # find empty columns
    my %empty;
    for my $c (@cols) {
	if (defined $tree) {
	    $empty{$c} = 1;
	    for my $species (@node_id) {
		$empty{$c} = $empty{$c} && !defined ($species_col2seq{$species}->[$c]);
	    }
	} else {
	    $empty{$c} = !defined ($col2seq[$c]);
	}
    }

    # eliminate empty columns (but ensure alignment is non-empty)
    @cols = map ($empty{$_} ? () : $_, @cols);
    if (@cols == 0) {
	@cols = (0);
	$col2ss[0] = $gapchar;  # "-"
    }

    # create & populate Stockholm object
    my $stock = Stockholm->new;
    $stock->gc_SS_cons (join ("", map ($col2ss[$_], @cols)));

    if (defined $tree) {
	for my $species (@node_id) {
	    my $col2seq = $species_col2seq{$species};
	    $col2seq = [] unless defined $col2seq;
	    $stock->add_row ($species, join ("", map (defined() ? $_ : "-", map ($col2seq->[$_], @cols))));
	}

	$stock->add_gf ("NH", $tree->to_string);
    } else {
	$stock->add_row ("seq", join ("", map (defined() ? $_ : "-", map ($col2seq[$_], @cols))));
    }

    # return
    return $stock;
}

# make_tkfst_string method
sub make_tkfst_string {
    my ($self) = @_;
    my $tkfst = "";
    $self->traverse
	(sub { # pre
	    my ($node) = @_;
	    $tkfst .= $node->state;
	    $tkfst .= "-" if @{$node->child} == 1;
	    $tkfst .= "(" if @{$node->child} > 1;
	 },
	 sub { # mid
	     my ($node, $n_child) = @_;
	     $tkfst .= ")" if $n_child < @{$node->child} - 1;
	     $tkfst .= "(" if $n_child < @{$node->child} - 2;
	 });
    return $tkfst;
}

# Method to store the set of mutations at a node.
# Sets $liveref->{$node} to a reference to an array, every element of which is a reference to a 2-element array: [$mutation_rate, $coderef_implementing_mutation].
# The $coderef's are closures on the entire TKFST doubly-linked tree structure, and they update themselves by callback to this method. (Neat, huh)
# If $node has no mutations, then this subroutine deletes the entry $liveref->{$node}, so that keys(%$liveref) can always be used to get a set of "live" nodes.
sub set_mutations {
    my ($node, $params, $log, $liveref) = @_;
    my $state = $node->state;

    my $update_live = sub {
	my ($node) = @_;
	$node->set_mutations ($params, $log, $liveref);
    };

    my @mutation;
    if (is_live ($state)) {

	# detect the special case that this is the S state at the root (which occurs when $s_root==1)
	my $is_sroot = $state eq $stemState && !defined ($node->parent);

	# deletions
	if ($state ne $loopState) {  # L's can never be directly deleted; all other states can
	    if (length ($state) == 1 && !$is_sroot) {  # includes S states, but not S at the root
		push @mutation, [$params->mu_loop,
				 sub {
				     &$log ("Loop delete", $node, @_);
				     if ($state eq $stemState) {
					 # remove substructure nodes from live set
					 $node->child->[0]->traverse
					     (sub {
						 my ($descendant) = @_;
						 delete $liveref->{$descendant};
						 $descendant->state ($deadState);
					      });
				     }
				     delete $liveref->{$node};  # remove node from live set
				     $node->state ($deadState);
				 }];  # deletion of a base or substructure in a loop sequence

	    } elsif (length ($state) == 2) {  # basepair
		push @mutation, [$params->mu_stem,
				 sub {
				     &$log ("Stem delete", $node, @_);
				     delete $liveref->{$node};  # remove node from live set
				     $node->state ($deadState);
				 }];  # deletion of a basepair in a stem sequence
	    }
	}

	# substitutions
	if ($state ne $loopState && $state ne $stemState) {  # base or basepair
	    my $mutref = $params->subrate->{$state};
	    for my $destState (keys %$mutref) {
		push @mutation, [$mutref->{$destState},
				 sub {
				     &$log ("Subst to $destState", $node, @_);
				     $node->state ($destState);
				     # update live mutation set
				     &$update_live ($node);  # update substitutions
				 }];  # substitution in a loop or stem sequence
	    }
	}

	# insertions
	if (length ($state) == 1) {  # includes L and S states, as well as bases
	    # handle insertions in the parent loop sequence
	    if (!$is_sroot) {  # exclude S's at the root, since they have no parent loop sequence
		push @mutation, [$params->lambda_loop,
				 sub {
				     &$log ("Loop insert", $node, @_);
				     my $next_node = pop @{$node->child};
				     my $inserted_node = sample_loop_node ($params, $node);
				     if (defined $next_node) {
					 $next_node->parent ($inserted_node);
					 push @{$inserted_node->child}, $next_node;
				     }
				     # update live mutation set
				     $inserted_node->traverse ($update_live);
				     &$update_live ($node);  # update insertion links
				 }];  # insertion of a base or substructure in a loop sequence
	    }

	    # handle insertions after the "immortal link" of stem sequences
	    # (actually "immortal" is a bit of a misnomer, because it can be deleted by the parent sequence)
	    if ($state eq $stemState) {
		push @mutation, [$params->lambda_stem,
				 sub {
				     &$log ("Start-of-stem insert", $node, @_);
				     my $next_in_loop = @{$node->child} > 1 ? pop @{$node->child} : undef;
				     my $next_in_stem = pop @{$node->child};
				     my $inserted_node = sample_stem_node ($params, $node);
				     $next_in_stem->parent ($inserted_node);
				     push @{$inserted_node->child}, $next_in_stem;
				     if (defined $next_in_loop) {
					 push @{$node->child}, $next_in_loop;
				     }
				     # update live mutation set
				     $inserted_node->traverse ($update_live);
				     &$update_live ($node);  # update insertion links
				 }];  # insertion of a basepair at the beginning of a stem sequence
	    }

	} elsif (length ($state) == 2) {  # basepair states
	    # handle insertions after "mortal links" of stem sequences
	    push @mutation, [$params->lambda_stem,
			     sub {
				 &$log ("Stem insert", $node, @_);
				 my $next_node = pop @{$node->child};
				 my $inserted_node = sample_stem_node ($params, $node);
				 if (defined $next_node) {
				     $next_node->parent ($inserted_node);
				     push @{$inserted_node->child}, $next_node;
				 }
				 # update live mutation set
				 $inserted_node->traverse ($update_live);
				 &$update_live ($node);  # update insertion links
			     }];  # insertion of a basepair in a stem sequence
	}
    }

    # update %{$liveref}
    if (@mutation) {
	$liveref->{$node} = \@mutation;
    } elsif (exists $liveref->{$node}) {
	delete $liveref->{$node};
    }
}

# method to evolve for a given branch length, using Gillespie's algorithm.
# returns an event log summary: an array of [TAG,SUMMARY] tuples that can be added to a Stockholm file as #=GF TAG SUMMARY
sub evolve {
    my ($self, $params, $time, $history_callback) = @_;

    # declare $time_left
    my $time_left;

    # logging function
    my %log;
    my $events = 0;
    my $log = sub {
	my ($message, $node, @args) = @_;
	++$log{$message}->{$node->state};
	++$events;
	my $time_now = $time - $time_left;
	# description of event
	print STDERR " $message at node state '", $node->state, "', time $time_now" if $verbose;
	# call to $history_callback
	if (defined $history_callback) {
	    &$history_callback ($time_now, $message, $node, @args);
	}
	# number of events
	warn "Logged $events events" if $events % 1000 == 0;
    };

    # for speed, build a map from "live" nodes to their possible mutations
    # (we could traverse the entire tree every step of the main loop, but scanning all the dead nodes takes time)
    my %live;
    my $update_live = sub {
	my ($node) = @_;
	$node->set_mutations ($params, $log, \%live);
    };
    $self->traverse ($update_live);

    # initial call to $history_callback
    if (defined $history_callback) {
	&$history_callback (0);
    }

    # main loop
    for ($time_left = $time; $time_left > 0; ) {

	# uncomment the following line to rebuild %live at every iteration
	# (slow, but useful for debugging)
#	$self->traverse ($update_live);

	# get total rate; also count number of live nodes & possible mutations
	my $totalrate = 0;
	my @nodes = keys %live;
	my $nmut = 0;
	for my $node (@nodes) {
	    my $node_is_live = 0;
	    for my $mut (@{$live{$node}}) {
		$totalrate += $mut->[0];
		++$nmut;
	    }
	}

	# debug/status message
	print STDERR "  (", @nodes+0, " live nodes; $nmut possible mutations; total mutation rate = $totalrate; time left = $time_left)\n" if $verbose;

	# sample wait time to next event
	my $wait_time = -log(randprob()) / $totalrate;
	last if $wait_time > $time_left;
	$time_left -= $wait_time;

	# sample nature of next event
	my $event_type = randprob() * $totalrate;
	my $event = sub { warn "Warning: undefined mutation event!" };
      NODELOOP: for my $node (@nodes) {
	    for my $mut (@{$live{$node}}) {
		if (($event_type -= $mut->[0]) <= 0) {
		    $event = $mut->[1];
		    last NODELOOP;
		}
	    }
	}

	# do next event
	&$event();

	# print Stockholm alignment (debug)
#	warn $self->make_stockholm->to_string;
    }
    print STDERR "\n" if $verbose;

    # final call to $history_callback
    if (defined $history_callback) {
	&$history_callback ($time);
    }

    # print event summary
    my @summary;
    push @summary, "Length $time.";
    if ($events) {
	for my $message (sort keys %log) {
	    push @summary,
	    "  $message: ",
	    join (", ",
		  map ("$_($log{$message}->{$_})",
		       sort keys %{$log{$message}})),
	    ".";
	}
    } else {
	push @summary, "  No events.";
    }

    # return summary
    return join ("", @summary);
}

# method to evolve on an entire tree
sub evolve_tree {
    my ($self, $params, $tree) = @_;

    # track event summaries
    my @tag_summary;

    # get "node IDs" for all nodes in Newick tree
    my @node_id = get_node_id ($tree);

    # evolve
    for (my $tree_node = 0; $tree_node < $tree->nodes; ++$tree_node)
    {
	my $nodeName = $node_id[$tree_node];
	# evolve TKFST along a branch, unless this is the root
	if ($tree_node == 0) {
	    push @tag_summary, ['ROOT', $self->make_tkfst_string];
	} else {
	    # log
	    my $parentName = $node_id[$tree->parent->[$tree_node]];
	    print STDERR "Evolving sequence $parentName into sequence $nodeName" if $verbose;

	    # get parent TKFST
	    $self->restore_state ($parentName);

	    # open history file & create callback sub
	    my $historyCallback;
	    local *HISTORY;
	    if ($history) {
		open HISTORY, ">$nodeName$historySuffix";
		$historyCallback = sub {
		    my ($time, $message, $node, @args) = @_;
		    # quick & dirty way to get sequence & structure strings, using make_stockholm method
		    $self->save_state ($nodeName);
		    my $tmpStock = $self->make_stockholm ($tree);
		    my $ss = $tmpStock->gc_SS_cons;
		    my $seq = $tmpStock->seqdata->{$nodeName};  # current state is stored under parentName
		    # copy gaps from $seq to $ss, then remove all gaps... ugh
		    for (my $pos = 0; $pos < length($seq); ++$pos) {
			substr($ss,$pos,1) = $gapchar if substr($seq,$pos,1) eq $gapchar;
		    }
		    $ss =~ s/$gapchar//g;
		    $seq =~ s/$gapchar//g;
		    warn $tmpStock->to_string, "\$nodeName = $nodeName\n\$parentName = $parentName\n\$ss  = '$ss'\n\$seq = '$seq'\nlength(\$ss) != length(\$seq)" if length($ss) != length($seq);
		    # make TKFST string
		    my $tkfst_string = $self->copy_live->make_tkfst_string;
		    # print log message
		    print HISTORY "$time $tkfst_string $seq $ss\n";
		};
	    }

	    # evolve
	    my $branch_length = $tree->branch_length->[$tree_node];
	    if (defined $branch_length) {
		my $summary = $self->evolve ($params, $branch_length, $historyCallback);
		push @tag_summary,
		['BRANCH', "$summary  ($parentName -> $nodeName)"],
		[$tree->children($tree_node) > 0 ? 'NODE' : 'LEAF', $self->make_tkfst_string];
	    }

	    # close history file
	    close HISTORY if $history;
	}
	
	# save state for children (and for make_stockholm)
	$self->save_state ($nodeName);
    }

    # return summaries
    return @tag_summary;
}

# get node IDs for a tree
sub get_node_id {
    my ($tree) = @_;
    my @id = @{$tree->node_name};
    return map (defined($id[$_]) && length($id[$_]) ? $id[$_] : ($_ == 0 ? "root" : "node$_"), 0..@id-1);
}


1;
