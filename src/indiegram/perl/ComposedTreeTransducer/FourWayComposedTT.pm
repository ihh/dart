#!/usr/bin/perl -w

package ComposedTreeTransducer::FourWayComposedTT;

use strict;
use Carp;
use Util;

use ComposedTreeTransducer;
use TreeTransducer::SingletTT;
use TreeTransducer::BranchTT;

our @ISA = qw /ComposedTreeTransducer/;
use vars '@ISA';

# To do (ordered high to low priority):
#  - Write _checkValid(), particularly something which checks that every dest state in the TM
#    is present as a source state as well.


# tree
my %parent = (
    2 => 1,
    3 => 1,
    4 => 1
    );
my %descendants = (
    1 => [qw/2 3 4/],
    2 => [qw//],
    3 => [qw//],
    4 => [qw//]
    );

# nodes
my $nodes = 4;
my @nodes = 1..$nodes;
my @leaves = 2..$nodes;

# branch lengths (keyed by leaf node)
my %branchLen = (
		 2 => "t",
		 3 => "u",
		 4 => "v"
		);

sub new {
  my ($class, $singletfile, $branchfile, $directives) = @_;

  my $self = ComposedTreeTransducer->new ($directives);

  bless $self, ref ($class) || $class;

  $self->{'singlet'} = "";
  $self->{'branch'} = "";
  $self->{'reducedStateGraph'} = {};
  $self->{'reducedTM'} = {};
  $self->{'fullyReducedStateGraph'} = {};
  $self->{'fullyReducedTM'} = {};
  $self->{'fullyReducedBifurcMap'} = {};

  $self->singlet (TreeTransducer::SingletTT->new ($singletfile));
  $self->branch (TreeTransducer::BranchTT->new ($branchfile));

  $self->_initialize();

  return $self;
}

sub _initialize {
  my ($self) = @_;  

  # store the "overall" composed Start state
  my $tmp = {};
  $tmp->{1} = $self->singlet->startState();
  foreach my $leaf (@leaves) {
    $tmp->{$leaf} = $self->branch->startState();
  }
  $self->startState ($tmp);

  # store the composed End state
  $tmp = {};
  $tmp->{1} = $self->singlet->endState();
  foreach my $leaf (@leaves) {
    $tmp->{$leaf} = $self->branch->endState();
  }
  $self->endState ($tmp);


  # Create and store the state graph.
  # The state graph holds information about accessible states for all states,
  # including bifurcations.  The list of accessible states for bifurcation states
  # is their child states.
  #   $self->stateGraph->{$srcKey} is an array ref of dest composed state keys
  # The all-End state is present in the state graph.
  $self->stateGraph ($self->buildStateGraph ($self->startState));

  # Create and store the "bifurcation transition matrix."
  #   $self->bifurcMap->{$srcKey}->{'left'}->{$childKey} is a string giving
  #   the corresponding transition probability
  # Bifurcations are deterministic here, so there is only 
  # one childKey and all probabilities are 1.
  $self->bifurcMap ($self->buildBifurcMap ($self->stateGraph));

  # Remove bifurcation states which have 2 End children.
  # Update the stateGraph and bifurcMap accordingly.
  unless ($self->directives_noprune) {
    $self->pruneBifurcations();
  }

  # Create and store the transition matrix.
  # The transition matrix holds information about transition probabilities 
  # between all valid source and destination states, excepting bifurcations.
  #   $self->tm->{$srcKey}->{$destKey} is a string giving
  #   the corresponding transition probability
  # NB: Bifurcations are NOT present as source states in
  # the transition matrix, and neither is the End state.
  $self->tm ($self->buildTM ($self->stateGraph));

  # Create and store the reduced state graph.
  # The reduced state graph is derived from the state graph by eliminating all
  # Null states ("windback" states), excepting only those which are children
  # of bifurcation states.
  $self->reducedStateGraph ($self->buildReducedStateGraph());

  # Create and store the reduced transition matrix.
  # (TM with all Null states removed, excepting those which are children of
  # bifurcation states).
  $self->reducedTM ($self->buildTM ($self->reducedStateGraph));

  # Create and store the fully reduced state graph.
  # The fully reduced state graph is derived from the state graph by eliminating all
  # Null states ("windback" states), including those which are children
  # of bifurcation states.
  $self->fullyReducedStateGraph ($self->buildFullyReducedStateGraph());

  # Create and store the fully reduced "bifurcation transition matrix".
  # Fully null state elimination generically gives rise to non-deterministic
  # bifurcations.
  #   $self->fullyReducedBifurcMap->{$srcKey}->{'left'}->{$childKey} is a string giving
  #   the corresponding transition probability
  $self->fullyReducedBifurcMap ($self->buildFullyReducedBifurcMap());

  # Create and store the fully reduced transition matrix
  # (TM with all Null states removed).
  $self->fullyReducedTM ($self->buildTM ($self->fullyReducedStateGraph));

}

sub _checkValid {
  my ($self) = @_;

  # check that all branch bifurc states have absorb profiles
  # matching a singlet bifurc emit profile
  foreach my $s (@{$self->branch->states()}) {
    if ($self->branch->isBifurcMatch ($s)) {
      my $ok = 0;
      foreach my $b (@{$self->singlet->states()}) {
	if ($self->singlet->isBifurcInsert ($b)) {
	  if ($self->areCompatBifurc ($b, $s)) { $ok = 1; }
	}
      }
      if (!$ok) { croak ("Warning: No singlet BifurcInsert state has an emit profile matching the absorb profile of (branch) state '$s'."); }
    }
  }

  # checks on transition matrix and state graph

  # checks on reduced tm and reduced state graph

  # check that probability is conserved by the transition matrix

}

################################################
# Singlet/branch emission compatibility checks #
################################################

# Do two state types of type Insert and Match have compatible
# emit and absorb profiles?
sub areCompatEmit {
  my ($self, $insertState, $matchState) = @_;

  if ($self->singlet->isInsert ($insertState) && $self->branch->isMatch ($matchState)) {
    if ($self->singlet->emitProfile ($insertState) eq $self->branch->absorbProfile ($matchState)) {
      return 1;
    }
  }
  return 0;
}

# Do two states of type Bifurc have compatible emit and absorb profiles?
sub areCompatBifurc {
  my ($self, $insertBifurcState, $matchBifurcState) = @_;

  if ($self->singlet->isBifurcInsert ($insertBifurcState) && $self->branch->isBifurcMatch ($matchBifurcState)) {
    if ($self->singlet->emitProfile ($insertBifurcState) eq $self->branch->absorbProfile ($matchBifurcState)) {
      return 1;
    }
  }
  return 0;
}

###################################
# Rules of transducer composition #
###################################

# Is a composed state a Start state?
sub isStart {
  my ($self, $cState) = @_;


  # 2 cases
  # case 1: bifurcation at the root, passed to some number of leaf nodes
  # (meaning all leaves at either Start or End)
  if ($self->singlet->isStart ($cState->{1})) {
    my $is = 1;
    foreach my $leaf (@leaves) {
      unless ($self->branch->isStart ($cState->{$leaf}) || $self->branch->isEnd ($cState->{$leaf})) { $is = 0; }
    }
    if ($is) { return 1; }
  }

  # case 2: bifurcation at a leaf
  # (meaning one leaf in Start state; all other nodes in End state)
  if ($self->singlet->isEnd ($cState->{1})) {
    my $startCnt = 0;
    my $endCnt = 0;
    foreach my $leaf (@leaves) {
      if ($self->branch->isStart ($cState->{$leaf})) { $startCnt += 1; }
      elsif ($self->branch->isEnd ($cState->{$leaf})) { $endCnt += 1; }
    }
    if (($startCnt eq 1) && ($endCnt eq (scalar(@leaves)-1))) { return 1; }
  }

  return 0;

}

# Is a composed state the (unique) composed End state?
sub isEnd {
  my ($self, $cState) = @_;
  my $is = 1;
  unless ($self->singlet->isEnd ($cState->{1})) { $is = 0; }
  foreach my $leaf (@leaves) {
    unless ($self->branch->isEnd ($cState->{$leaf})) { $is = 0; }
  }
  return $is;
}

# Given a composed state, is it an emit state?
sub isEmit {
  my ($self, $cState) = @_;

  my $activeNode = $self->activeNode ($cState);

  # 2 cases:
  # case 1: root in Insert state; all other nodes in End state
  if ($activeNode == 1) {
    my $is = 1;
    foreach my $leaf (@leaves) {
      unless ($self->branch->isEnd ($cState->{$leaf})) { $is = 0; }
    }
    return ($is && $self->singlet->isInsert ($cState->{1}));
  }
  # case 2: composed Emit state
  elsif ($activeNode > 1 && $activeNode <= $nodes) {
    # 2 possible cases for composed emit states
    #  NB: must be considered in this order!
    # case 1: insertion to single leaf
    if ($self->branch->isInsert ($cState->{$activeNode})) {
      return 1;
    }
    # case 2: insertion at root, passed to leaves
    #  root in Insert state, all leaves in Match or End state
    elsif ($self->singlet->isInsert ($cState->{1})) {
      foreach my $leaf (@leaves) {
	unless ($self->branch->isMatch ($cState->{$leaf}) || $self->branch->isEnd ($cState->{$leaf})) { return 0; }
      }
      return 1;
    }
  }

  return 0;
}

# Given a composed state, is it a null state?
sub isNull {
  my ($self, $cState) = @_;
  return (!($self->isStart ($cState) || $self->isEnd ($cState) || $self->isEmit ($cState) || $self->isBifurc ($cState)));
}

# Given a composed state, is it a bifurcation state?
# Note that this covers the case of the bifurcation originating at the
# root as well as at a leaf node.
sub isBifurc {
  my ($self, $cState) = @_;

  # 2 cases
  # case 1: bifurcation originated at the root
  my $is = 1;
  unless ($self->singlet->isBifurc ($cState->{1})) { $is = 0; }
  foreach my $leaf (@leaves) {
    unless ($self->branch->isBifurc ($cState->{$leaf}) || $self->branch->isEnd ($cState->{$leaf})) { $is = 0; }
  }
  if ($is) { return 1; }

  # case 2: bifurcation originated at a leaf
  unless ($self->singlet->isBifurc ($cState->{1})) {
    my $bifurcCnt = 0;
    foreach my $leaf (@leaves) {
      if ($self->branch->isBifurc ($cState->{$leaf})) { $bifurcCnt += 1; }
    }
    if ($bifurcCnt eq 1) { return 1; }
    elsif ($bifurcCnt > 1) { croak ("state '", $self->cStateKey ($cState), " isn't a valid composed state.\n"); }
  }

  return 0;
}

# Returns a string {Start, End, Null, Emit, Bifurc}.
sub stateType {
  my ($self, $cStateKey) = @_;

  my $cState = $self->cState ($cStateKey);
  unless (defined $cState) { die; }
  if ($self->isStart ($cState)) { return 'Start'; }
  elsif ($self->isEnd ($cState)) { return 'End'; }
  elsif ($self->isNull ($cState)) { return 'Null'; }
  elsif ($self->isEmit ($cState)) { return 'Emit'; }
  elsif ($self->isBifurc ($cState)) { return 'Bifurc'; }
  else { croak "Unreachable.\n"; }
}

# State comparison for sorting.
# Input is two composed state keys.
sub stateCmp {
  my ($self, $aKey, $bKey) = @_;
  my $a = $self->cState ($aKey);
  my $b = $self->cState ($bKey);

  # root dominant
  my $score = ($self->singlet->stateSorting->{$a->{1}} <=> $self->singlet->stateSorting->{$b->{1}}) * 1000;

  # then weight in favor of postorder updates
  foreach my $leaf (@leaves) {
    $score += ($self->branch->stateSorting->{$a->{$leaf}}*2 <=> $self->branch->stateSorting->{$b->{$leaf}}*2) * $leaf;
  }

  return $score;
}


# Given a composed state, get the active node.
# Return a dummy value -1 if all-End state.
sub activeNode {
  my ($self, $cState) = @_;

  # do a postorder traversal of tree, looking for non-Wait and non-End nodes

  # leaves
  foreach my $n (reverse (@leaves)) { # postorder
    my $s = $cState->{$n};
    unless ($self->branch->isWait ($s) || $self->branch->isEnd ($s)) {
      return $n;
    }
  }
  # then root
  unless ($self->singlet->isEnd ($cState->{1})) {
    return 1;
  }
  # if we reach here, then it must be the all-End state
  if ($self->isEnd ($cState)) {
    return -1;
  } else {
    croak ("Unreachable code: state ", $self->cStateKey ($cState), " has no active node.\n");
  }
}

# Given a composed state of type Emit, get the emitting node
# (this is the highest-numbered node of type Insert).
# Returns the dummy value -1 if state not of type Emit.
sub emitter {
  my ($self, $cState) = @_;

  # return dummy value if not Emit state
  unless ($self->isEmit ($cState)) {
    return -1;
  }

  # Note that this is effectively equivalent to performing a postorder traversal
  # and returning the first insert state found.
  my $node = -1;
  # set $node to root if in an emit state
  if ($self->singlet->isInsert ($cState->{1})) {
    $node = 1;
  }
  foreach my $leaf (@leaves) {
    if ($self->branch->isInsert ($cState->{$leaf})) {
      $node = $leaf;
    }
  }

  return $node;
}

# Given a composed state of type Emit, get the set of nodes in emit states.
# Returns an empty array if state not of type Emit.
sub emitSet {
  my ($self, $cState) = @_;

  # if not an Emit state, return an empty array
  unless ($self->isEmit ($cState)) { 
    return qw();
  }
    
  my $emitter = $self->emitter ($cState);
  my @emitSet;
  if ($emitter eq 1) {
    push (@emitSet, $emitter) if ($self->singlet->isInsert ($cState->{$emitter}));
    foreach my $n (@leaves) {
      push (@emitSet, $n) unless ($self->branch->isEnd ($cState->{$n}));
    }
  } elsif ($emitter > 1 && $emitter <= $nodes) {
    push (@emitSet, $emitter) if ($self->branch->isInsert ($cState->{$emitter}));
  }
  
  # sanity check: there had better be a non-empty emitset
  unless (scalar (@emitSet)) {
    croak ($self->cStateKey ($cState), " has an empty emit set.\n");
  }

  return @emitSet;
}

# Given source and destination composed states, return 
# the set of nodes which have changed state during the transition.
# This doesn't really do robust checking that that src and dest states are
# valid or that the transition itself is valid; hopefully the upstream
# code will only call it on valid transitions.
# This works for transitions defined on both the un-reduced and reduced state graphs/TMs.
sub mutableNodes {
  my ($self, $src, $dest) = @_;

  my @mutable;

  # NB: The case of a Null dest state occurs only for the un-reduced TM.
  # (the mutable nodes are just the active node of the source state plus its non-End descendants).

  # Transitions from bifurcations must be handled first!
  if ($self->isBifurc ($src)) {
    # case 1: bifurcation originated at root
    if ($self->singlet->isBifurc ($src->{1})) {
      push (@mutable, 1);
      foreach my $leaf (@leaves) {
	push (@mutable, $leaf) unless ($self->branch->isEnd ($src->{$leaf}));
      }
    }
    # case 2: bifurcation originated at leaf
    else {
      foreach my $leaf (@leaves) {
	if ($self->branch->isBifurc ($src->{$leaf})) {
	  push (@mutable, $leaf);
	}
      }
      if (scalar (@mutable) > 1) { croak "Too many bifurc states!  Not a valid composed state.\n"; }
    }
  }
  elsif ($self->isNull ($dest)) {
    my $activeNode = $self->activeNode ($src);
    push (@mutable, $activeNode);
    foreach my $n (@{$descendants{$activeNode}}) {
      unless ($self->branch->isEnd ($src->{$n})) {
	push (@mutable, $n);
      }
    }
  }
  elsif ($self->isEmit ($dest)) {
    # case 1: emission at root
    my $emitter = $self->emitter ($dest);
    if ($emitter eq 1) {
      push (@mutable, $self->emitSet ($dest));
    }
    # case 2: emission at a leaf
    # store all updated nodes (meaning all nodes encountered before $emitter in a
    # post-order traversal which aren't in Wait or End states)
    elsif ($emitter > 1 && $emitter <= $nodes) {
      foreach my $leaf ($emitter..$nodes) {
	my $s = $src->{$leaf};
	unless ($self->branch->isWait ($s) || $self->branch->isEnd ($s)) {
	  push (@mutable, $leaf);
	}
      }
    } else { croak "Unreachable code\n"; }
  }
  elsif ($self->isBifurc ($dest)) {
    # case 1: bifurcation originated at root
    if ($self->singlet->isBifurc ($dest->{1})) {
      push (@mutable, 1);
      foreach my $leaf (@leaves) {
	push (@mutable, $leaf) unless ($self->branch->isEnd ($dest->{$leaf}));
      }
    }
    # case 2: bifurcation originated at leaf; only that leaf changed state
    else {
      foreach my $leaf (@leaves) {
	if ($self->branch->isBifurc ($dest->{$leaf})) {
	  push (@mutable, $leaf);
	}
      }
    }
  }
  elsif ($self->isStart ($dest)) {
    # case 1: src state has a bifurcation at the root
    # meaning all dest states are in Start or End states
    if ($self->singlet->isStart ($dest->{1})) {
      push (@mutable, 1);
      foreach my $leaf (@leaves) {
	push (@mutable, $leaf) unless ($self->branch->isEnd ($dest->{$leaf}));
      }
    }
    # case 2: src state has a BifurcInsert at a leaf,
    # and this is the left or right composed child (meaning 
    # all nodes except the bifurcating one are in End states)
    elsif ($self->singlet->isEnd ($dest->{1})) {
      foreach my $leaf (@leaves) {
	if ($self->branch->isStart ($dest->{$leaf})) {
	  push (@mutable, $leaf);
	}
      }
    }
    else { croak "Unreachable code\n"; }
  }
  elsif ($self->isEnd ($dest)) {
    # case 1: transition to End originated at root
    unless ($self->singlet->isEnd ($src->{1})) {
      push (@mutable, 1);
      foreach my $leaf (@leaves) {
	push (@mutable, $leaf) unless ($self->branch->isEnd ($src->{$leaf}));
      }
    }
    # case 2: transition to End originated at leaf
    else {
      foreach my $leaf (@leaves) {
	unless ($self->branch->isEnd ($src->{$leaf})) {
	  push (@mutable, $leaf);
	}
      }
    }
  }
  else {
    croak "Unreachable code.\n";
  }

  # sanity check
  unless (scalar (@mutable)) { croak ("Error, no mutable nodes for transition", $self->cStateKey ($src), " -> ", $self->cStateKey ($dest), "\n"); }

  return @mutable;
}


# Given a composed state, return a list of references to
# the allowed destination null composed states.
# The state of a single node is updated.
sub allowedNullTrans {
  my ($self, $cState) = @_;

  my $activeNode = $self->activeNode ($cState);

  # results: allowed destination composed states
  my @destComposedStates;

  # root
  if ($activeNode eq 1) {
    foreach my $destState (@{$self->singlet->destStates ($cState->{$activeNode})}) {
      
      # if null transition, store
      if ($self->singlet->isStart ($destState)) {
	my %destComposedState = %$cState;
	$destComposedState{$activeNode} = $destState;
	push (@destComposedStates, \%destComposedState);
      }
      
    }
  }
  # leaves
  elsif ($activeNode > 1 && $activeNode <= $nodes) {
    foreach my $destState (@{$self->branch->destStates ($cState->{$activeNode})}) {
      
      # if null transition, store
      if ($self->branch->isStart ($destState) || $self->branch->isWait ($destState)) {
	my %destComposedState = %$cState;
	$destComposedState{$activeNode} = $destState;
	push (@destComposedStates, \%destComposedState);
      }
      
    }
  }
  
  return @destComposedStates;
}

# Given a composed state, return a list of references to
# the allowed destination emit composed states.
# (The list of allowed destination states includes only
# those where the emission has been "fully processed,",
# meaning propagated down the tree.)
sub allowedEmitTrans {
  my ($self, $cState) = @_;

  # allowed destination composed states
    my @destComposedStates;
    my $activeNode = $self->activeNode ($cState);

    if ($activeNode eq 1) {

	# loop over all destination states for the root node
	foreach my $rootDestState (@{$self->singlet->destStates ($cState->{$activeNode})}) {  # $rootDestState is the active node destination state

	    # we're only interested in transitions to emit states
	    unless ($self->singlet->isInsert ($rootDestState)) { next; }

	    # sanity check that leaves are in Wait or End states
	    foreach my $leaf (@leaves) {
	      unless ($self->branch->isWait ($cState->{$leaf}) || $self->branch->isEnd ($cState->{$leaf}))
		{ croak ("Error: leaf not in Wait or End state in '", $self->cStateKey ($cState), "'\n"); }
	    }

	    # update the root node state	    
	    my %destComposedState = %$cState;
	    $destComposedState{$activeNode} = $rootDestState;

	    # An ugly hack:
	    # The below nested loops won't fire properly if one of the outer loops
	    # (leaf 2 or 3) is over an empty array.  The End state has no destStates,
	    # so states like S_e_e_WS are problematic.  We "fix" this by 
	    # pushing on the End state if the array will be empty.
	    my @leaf2DestStates;
	    push (@leaf2DestStates, $self->branch->isEnd ($cState->{2}) ? $cState->{2} : @{$self->branch->destStates ($cState->{2})});
	    my @leaf3DestStates;
	    push (@leaf3DestStates, $self->branch->isEnd ($cState->{3}) ? $cState->{3} : @{$self->branch->destStates ($cState->{3})});
	    my @leaf4DestStates;
	    push (@leaf4DestStates, $self->branch->isEnd ($cState->{4}) ? $cState->{4} : @{$self->branch->destStates ($cState->{4})});
	    
	    # loop over all combinations of leaf states which have appropriate
	    # absorption profiles
	    foreach my $leaf2DestState (@leaf2DestStates) {
	      foreach my $leaf3DestState (@leaf3DestStates) {
		foreach my $leaf4DestState (@leaf4DestStates) {
		  if (($self->branch->isEnd ($cState->{2}) ? 1 : $self->areCompatEmit ($rootDestState, $leaf2DestState)) and
		      ($self->branch->isEnd ($cState->{3}) ? 1 : $self->areCompatEmit ($rootDestState, $leaf3DestState)) and
		      ($self->branch->isEnd ($cState->{4}) ? 1 : $self->areCompatEmit ($rootDestState, $leaf4DestState))) {
			    my %tmp = %destComposedState;
			    $tmp{2} = $leaf2DestState;
			    $tmp{3} = $leaf3DestState;
			    $tmp{4} = $leaf4DestState;
			    push (@destComposedStates, \%tmp);
			}
		    }
		}
	    }

	}

    }

    # if emitter is leaf node, then just need to update the state for that branch
    # (because no descendants)
    elsif ($activeNode > 1 && $activeNode <= $nodes) {

      # loop over destination states
	foreach my $leafDestState (@{$self->branch->destStates ($cState->{$activeNode})}) {

	  # we're only interested in transitions to emit states
	  if ($self->branch->isInsert ($leafDestState)) {
		# update the active (leaf) node
		my %destComposedState = %$cState;
		$destComposedState{$activeNode} = $leafDestState;

		# store the new state
		push (@destComposedStates, \%destComposedState);
	    }
	}
    }

    return @destComposedStates;
}

# Given a composed state, return a list of references to the allowed destination
# bifurcation composed states.  (The list of allowed destination states
# includes only those where the emission has been "fully processed,"
# meaning propagated down the tree.)
sub allowedBifurcTrans {
  my ($self, $cState) = @_;

  # allowed destination composed states
  my @destComposedStates;
  my $activeNode = $self->activeNode ($cState);

  if ($activeNode eq 1) {
    
    # loop over destination states for the root node
    foreach my $rootDestState (@{$self->singlet->destStates ($cState->{$activeNode})}) {  # $rootDestState is the active node destination state

      # we're only interested in bifurcations
      unless ($self->singlet->isBifurc ($rootDestState)) { next; }

      # sanity check that leaves are in Wait or End states
      foreach my $leaf (@leaves) {
	unless ($self->branch->isWait ($cState->{$leaf}) || $self->branch->isEnd ($cState->{$leaf}))
	  {  croak ("Error: leaf not in Wait or End state in '", $self->cStateKey ($cState), "'\n"); }
      }

      # update the root node state	    
      my %destComposedState = %$cState;
      $destComposedState{$activeNode} = $rootDestState;

      # An ugly hack:
      # The below nested loops won't fire properly if one of the outer loops
      # (leaf 2 or 3) is over an empty array.  The End state has no destStates,
      # so states like S_e_e_WS are problematic.  We "fix" this by 
      # pushing on the End state if the array will be empty.
      my @leaf2DestStates;
      push (@leaf2DestStates, $self->branch->isEnd ($cState->{2}) ? $cState->{2} : @{$self->branch->destStates ($cState->{2})});
      my @leaf3DestStates;
      push (@leaf3DestStates, $self->branch->isEnd ($cState->{3}) ? $cState->{3} : @{$self->branch->destStates ($cState->{3})});
      my @leaf4DestStates;
      push (@leaf4DestStates, $self->branch->isEnd ($cState->{4}) ? $cState->{4} : @{$self->branch->destStates ($cState->{4})});
	    
      # loop over all combinations of leaf states which have appropriate
      # absorption profiles
      foreach my $leaf2DestState (@leaf2DestStates) {
	foreach my $leaf3DestState (@leaf3DestStates) {
	  foreach my $leaf4DestState (@leaf4DestStates) {
	    if (($self->branch->isEnd ($cState->{2}) ? 1 : $self->areCompatBifurc ($rootDestState, $leaf2DestState)) and
		($self->branch->isEnd ($cState->{3}) ? 1 : $self->areCompatBifurc ($rootDestState, $leaf3DestState)) and
		($self->branch->isEnd ($cState->{4}) ? 1 : $self->areCompatBifurc ($rootDestState, $leaf4DestState))) {
	      my %tmp = %destComposedState;
	      $tmp{2} = $leaf2DestState;
	      $tmp{3} = $leaf3DestState;
	      $tmp{4} = $leaf4DestState;
	      push (@destComposedStates, \%tmp);
	    }
	  }
	}
      }

    }

  }

  # if active node is leaf node, then just need to update the state for that branch
  # (because no descendants)
  elsif ($activeNode > 1 && $activeNode <= $nodes) {
    
      # loop over destination states
	foreach my $leafDestState (@{$self->branch->destStates ($cState->{$activeNode})}) {

	  # we're only interested in BifurcInsert states
	  if ($self->branch->isBifurcInsert ($leafDestState)) {
		# update the active (leaf) node
		my %destComposedState = %$cState;
		$destComposedState{$activeNode} = $leafDestState;

		# leave other nodes unchanged
		# bifurcComposedChildren() will take care of "expanding" the bifurcation
		# and matching up nonterminals as appropriate

		# store the new state
		push (@destComposedStates, \%destComposedState);
	    }
	}

  }

  return @destComposedStates;
}

# Given a composed state of type bifurcation, return a 
# hash keyed by (left, center, right) whose values are
# references to the corresponding (composed) "child" states.
#   $children->{'left'} is the composed state ref for the left composed child
# Note that we only allow deterministic bifurcations.
# If the composed state is not of type bifurc, return an empty array.
sub bifurcComposedChildren {
  my ($self, $cState) = @_;

  # if not a bifurcation
  unless ($self->isBifurc ($cState)) { return qw(); }

  my $children;
  my %lchild;
  my %cchild;
  my %rchild;

  # case 1: bifurcation originating at root
  if ($self->singlet->isBifurc ($cState->{1})) {
    ($lchild{1}, $cchild{1}, $rchild{1}) = @{$self->singlet->destStates ($cState->{1})};
    foreach my $leaf (@leaves) {
      if ($self->branch->isEnd ($cState->{$leaf})) {
	($lchild{$leaf}, $cchild{$leaf}, $rchild{$leaf}) = ($self->branch->endState(), $self->branch->endState(), $self->branch->endState());
      }	else {
	($lchild{$leaf}, $cchild{$leaf}, $rchild{$leaf}) = @{$self->branch->destStates ($cState->{$leaf})};
      }
    }
  }

  else {
    # root
    # center branch of parse tree is "original"
    $cchild{1} = $cState->{1};
    # left and right branches: no nonterminals emitted, so go to End
    $lchild{1} = $rchild{1} = $self->singlet->endState();

    # leaves
    foreach my $leaf (@leaves) {
      # if there's a bifurcation here
      if ($self->branch->isBifurc ($cState->{$leaf})) {
	($lchild{$leaf}, $cchild{$leaf}, $rchild{$leaf}) = @{$self->branch->destStates ($cState->{$leaf})};
      }
      # else just like root
      else {
	$cchild{$leaf} = $cState->{$leaf};
	$lchild{$leaf} = $rchild{$leaf} = $self->singlet->endState();
      }
    }
  }

  if ($self->directives_debug) {
    print $self->cStateKey ($cState), " -> (", $self->cStateKey (\%lchild), ", ", $self->cStateKey (\%cchild), ", ", $self->cStateKey (\%rchild), ")\n";
  }

  $children->{'left'} = \%lchild;
  $children->{'center'} = \%cchild;
  $children->{'right'} = \%rchild;

  return $children;
}

# Given a composed state, is a transition to End allowed?
sub isAllowedEndTrans {
  my ($self, $cState) = @_;

  my $activeNode = $self->activeNode ($cState);

  # 2 cases:
  # case 1: transition to End initiated at root of tree
  if ($activeNode eq 1) {
    # look for a transition to End at root
    foreach my $rootState (@{$self->singlet->destStates ($cState->{1})}) {

      if ($self->singlet->isEnd ($rootState)) {
	# confirm that exist transitions to End from each non-End leaf state
	foreach my $leaf (@leaves) {
	  if ($self->branch->isEnd ($cState->{$leaf})) { next; }
	  my $allowed = 0;
	  foreach my $leafState (@{$self->branch->destStates ($cState->{$leaf})}) {
	    if ($self->branch->isEnd ($leafState)) {
	      $allowed = 1;
	      last;
	    }
	  }
	  if ($allowed eq 0) { return 0; }
	}
	return 1;
      }

    }
  }

  # case 2: bifurcation at a leaf node, meaning all nodes but one leaf in End state
  elsif ($activeNode > 1 && $activeNode <= $nodes) {
    # check that all nodes other than the active node are in End states
    unless ($self->singlet->isEnd ($cState->{1})) { return 0; }
    foreach my $leaf (@leaves) {
      if ($leaf eq $activeNode) { next; }
      unless ($self->branch->isEnd ($cState->{$leaf})) { return 0; }
    }
    # look for a transition to End at the active node
    foreach my $dest (@{$self->branch->destStates ($cState->{$activeNode})}) {
      if ($self->branch->isEnd ($dest)) { return 1; }
    }
  }

  return 0;
}

# Given 2 composed state keys, return the corresponding effective transition probability
# (transition probability on the reduced transition matrix built by buildReducedTM())
# as a formatted string.  If the transition is already defined in the full/original 
# transition matrix, then that transition probability is returned.
# This method defines all transition probabilities on the composed machine.
# Fancy stuff (bifurcations, etc.) is handled by TreeTransducer.
# Individual singlet/branch transitions are checked for validity, but
# the overall composed state transition is not.
sub effTransProb {
  my ($self, $srcKey, $destKey) = @_;

  my $src = $self->cState ($srcKey);
  my $dest = $self->cState ($destKey);

  # array of probability factors
  my @factors;

  # get the set of nodes which have changed state
  my @mutable = $self->mutableNodes ($src, $dest);
  foreach my $n (@mutable) {
    # root
    if ($n eq 1) { push (@factors, $self->singlet->effTransProb ($src->{$n}, $dest->{$n})); }
    # leaves
    else { push (@factors, $self->branch->effTransProb ($src->{$n}, $dest->{$n}, $branchLen{$n})); }
  }

  # return formatted string
#  return join (" * ", map ("(".$_.")", @factors));
  return join (" * ", @factors);
}

# Get the emit probability of a composed state.
# Returns a hash mapping each node in the emit set to its associated factor.
sub emitProbs {
  my ($self, $cState) = @_;
  unless ($self->isEmit ($cState)) { croak "'$cState' is not an Emit state.\n"; }

  my $factors = {};

  my @emitSet = $self->emitSet ($cState);
  foreach my $n (@emitSet) {
    my $s = $cState->{$n};
    if ($n eq 1) {
      $factors->{$n} = $self->singlet->emitDist->{$s} . "(" . $self->singlet->emitProfile ($s) . "_$n)";
    }
    elsif ($n > 1 && $n <= $nodes) {
      if ($self->branch->isEnd ($s)) { next; }
      elsif ($self->branch->isInsert ($s)) {
	$factors->{$n} = $self->branch->formatEmitDist ($self->branch->emitDist->{$s}, $n) . "(" . $self->branch->emitProfile ($s) . "_$n)"; }
      elsif ($self->branch->isMatch ($s)) {
	unless ($self->branch->emitProfile ($s) eq "") { # ignore delete states
	  $factors->{$n} = $self->branch->formatEmitDist ($self->branch->emitDist->{$s}, $n) . "(" . $self->branch->emitProfile ($s) . "_$n|" . $self->singlet->emitProfile ($cState->{1}) . "_1)";
	} }
      else { croak "Unreachable code.\n"; }
    }
    else { croak "Unreachable code.\n"; }
  }

  return $factors;
}

# Given a composed state, return a hash mapping emit nodes to their emit profiles.
sub emitProfile {
  my ($self, $cState) = @_;

  my $emitProfile = {};
  foreach my $n ($self->emitSet ($cState)) {
    if ($n eq 1) {
      $emitProfile->{$n} = $self->singlet->emitProfile ($cState->{$n});
    } else {
      $emitProfile->{$n} = $self->branch->emitProfile ($cState->{$n});
    }
  }

  return $emitProfile;
}

# Checks whether 2 composed states are equivalent.
# Helper for findRedundant().
sub areEquivalent {
  my ($self, $cKey1, $cKey2, $tm) = @_;

  my $cState1 = $self->cState ($cKey1);
  my $cState2 = $self->cState ($cKey2);

  my $t1 = $self->stateType ($cKey1);
  my $t2 = $self->stateType ($cKey2);

  unless ($t1 eq $t2) { return 0; }

  # check transitions and transition probs
  foreach my $destKey (keys %{$tm->{$cKey1}}) {
    my $found = 0;
    foreach my $candDestKey (keys %{$tm->{$cKey2}}) {
      if ($destKey eq $candDestKey and
	  $tm->{$cKey1}->{$destKey} eq $tm->{$cKey2}->{$candDestKey})
	{ $found = 1; }
    }
    unless ($found) { return 0; }
  }
  foreach my $destKey (keys %{$tm->{$cKey2}}) {
    my $found = 0;
    foreach my $candDestKey (keys %{$tm->{$cKey1}}) {
      if ($destKey eq $candDestKey and
	  $tm->{$cKey2}->{$destKey} eq $tm->{$cKey1}->{$candDestKey})
	{ $found = 1; }
    }
    unless ($found) { return 0; }
  }

  # now check emission profiles and probs for Emit states
  if ($t1 =~ /Emit/) {
    my $emitProb1 = $self->emitProbs ($cState1);
    my $emitProb2 = $self->emitProbs ($cState2);
    while (my ($n, $p) = each %$emitProb1) {
      unless (exists $emitProb2->{$n} and
	      $p eq $emitProb2->{$n})
	{ return 0; }
    }
    while (my ($n, $p) = each %$emitProb2) {
      unless (exists $emitProb1->{$n} and
	      $p eq $emitProb1->{$n})
	{ return 0; }
    }
  }

  return 1;
}


####################
# Graph algorithms #
####################

# Given a starting composed state %beginState (which will usually be $self->startState),
# recursively build the state graph:
# Get all unique composed states (eventually) accessible from the %beginState.
sub buildStateGraph {
  my ($self, $beginState) = @_;

  # results hash: map states to all accessible destination states
  #   $self->stateGraph->{$srcKey} is an array ref of dest composed state keys
  my $stateGraph = {};

  # bookkeeping hash: to prevent double-counting of states
  #  keyed by the composed state keys of accessible states
  #  an entry for a state key means that the state has been processed
  my $statesDone = {};
  $self->buildStateGraphHelper ($beginState, $statesDone, $stateGraph);

  # sanity check: bookkeeping hash and results hash should have same keys
  foreach my $x (keys %$statesDone) {
    unless (exists $stateGraph->{$x}) { 
      print ("Major error:\n");
      print "\$statesDone has '$x'; \$stateGraph doesn't.\n";
      print "\$stateGraph contains ";
      foreach my $s (keys %$stateGraph) {
	print " $s ";
      }
      print "\nThis could be because '$x' is an accessible state but has no destination states; check your grammar!";
      croak "Uh-oh!  Bookkeeping gone wrong.\n";
    }
  }

  return $stateGraph;
}

# Recursive helper: For each new destination state, recursively examine
# and store its child states.
# statesDone keeps tracks of states already examined (prevents infinite 
# recursion; remember, we're doing graph search!)
sub buildStateGraphHelper {
  my ($self, $beginState, $statesDone, $stateGraph) = @_;

  # if we've been here before, then return (prevent infinite recursion)
  if (defined $statesDone->{$self->cStateKey ($beginState)}) { return; }

  # update bookkeeping hash
  if ($self->directives_debug) { print "storing ", $self->cStateKey ($beginState), " in statesDone\n"; }
  $statesDone->{$self->cStateKey ($beginState)} = 1;

  # handle special case of End state separately
  #  (it has no destination states so the below loop won't catch it)
  if ($self->isEnd ($beginState)) {
    # store dummy array reference
    $stateGraph->{$self->cStateKey ($beginState)} = [qw //];
  }

  # get accessible destination composed states (Null, Emit and End)
  my @destStates; # array of composed dest state refs

  # if a bifurc state, dests are child states
  if ($self->isBifurc ($beginState)) {
    my $children = $self->bifurcComposedChildren ($beginState);
    foreach my $branch (keys %{$children}) { push (@destStates, $children->{$branch}); }
  }
  # otherwise dests are transitions to states of type Null, Emit, Bifurc and End
  else {
    push (@destStates, $self->allowedNullTrans ($beginState));
    push (@destStates, $self->allowedEmitTrans ($beginState));
    push (@destStates, $self->allowedBifurcTrans ($beginState));
    if ($self->isAllowedEndTrans ($beginState)) { push (@destStates, $self->endState); }
  }

  # loop over dest states
  foreach my $destState (@destStates) {

    # Store ref to dest state in the results hash
    # NB: Child states of bifurcations are stored as different possible destination states.
    if ($self->directives_debug) {
      print "storing ", $self->cStateKey ($destState), " in stateGraph->", $self->cStateKey ($beginState), "\n";
    }
    push (@{$stateGraph->{$self->cStateKey ($beginState)}}, $self->cStateKey ($destState));

    # call recursively
    $self->buildStateGraphHelper ($destState, $statesDone, $stateGraph);
  }

}



# Build the transition matrix corresponding to the passed state graph.
sub buildTM {
  my ($self, $stateGraph) = @_;

  my $tm = {};
  # now coerce this into the format we want (a la what is used for the singlet and branch trans matrices):
  #  ->{$src}->{$dest} is a transition probability string
  foreach my $srcKey (keys %{$stateGraph}) {
    # ignore bifurcations here
    if ($self->isBifurc ($self->cState ($srcKey))) { next; }

    foreach my $destKey (@{$stateGraph->{$srcKey}}) {
      $tm->{$srcKey}->{$destKey} = $self->effTransProb ($srcKey, $destKey);
    }
  }
 
  return $tm;
}

# Build the bifurcation map corresponding to the passed state graph.
sub buildBifurcMap {
  my ($self, $stateGraph) = @_;
  
  my $bifurcMap = {};
  foreach my $srcKey (keys %{$stateGraph}) {
    unless ($self->isBifurc ($self->cState ($srcKey))) { next; }

    my $children = $self->bifurcComposedChildren ($self->cState ($srcKey));
    foreach my $branch (keys %{$children}) {
      my $stateKey = $self->cStateKey ($children->{$branch});
      $bifurcMap->{$srcKey}->{$branch}->{$stateKey} = $self->effTransProb ($srcKey, $stateKey);
    }
  }

  return $bifurcMap;
}

# Identify bifurcations which have 2 End children.  Remove them from
# the tm and stateGraph.
sub pruneBifurcations {
  my ($self) = @_;

  # record which bifurcations are to be pruned
  my $prune = {};
  foreach my $srcKey (keys %{$self->bifurcMap}) {
    # loop over left, center, right children
    my @children;
    my $cntEnd = 0;
    foreach my $branch (keys %{$self->bifurcMap->{$srcKey}}) {
      my $childKey = (keys %{$self->bifurcMap->{$srcKey}->{$branch}})[0];
      if (!$self->isEnd ($self->cState ($childKey))) { push (@children, $childKey); }
      else { $cntEnd += 1; }
    }
    if ($cntEnd == 2) {
      unless (scalar (@children) == 1) { croak "Unreachable!\n"; }
      $prune->{$srcKey} = $children[0];
    }
  }

  # debugging
  if ($self->directives_debug) {
    print scalar (keys %$prune), " bifurc states to be pruned:\n";
    foreach my $p (keys %$prune) {
      print "  $p -> ", $prune->{$p}, "\n";
    }
  }

  # build the new (pruned) stateGraph
  my $prunedStateGraph = {};
  foreach my $srcKey (keys %{$self->stateGraph}) {
    # if a bifurc to be pruned, don't add it as a source state
    if (exists $prune->{$srcKey}) { next; }
    
    foreach my $destKey (@{$self->stateGraph->{$srcKey}}) {
      if (exists $prune->{$destKey}) {
	push (@{$prunedStateGraph->{$srcKey}}, $prune->{$destKey});
      } else {
	push (@{$prunedStateGraph->{$srcKey}}, $destKey);
      }
    }
  }

  # store the pruned stateGraph
  $self->stateGraph ($prunedStateGraph);

  # build and store the corresponding pruned bifurcMap
  $self->bifurcMap ($self->buildBifurcMap ($self->stateGraph));

}


# Returns a hash mapping Null state keys to array of refs to 
# all non-Null states to which they can transition either
#  1. directly
#  2. via other Null States
sub buildNullConnections {
  my ($self, $stateGraph) = @_;
  
  # hash mapping null state keys to non-null states linked by all-Null transitions
  my $nullConnections = {};
  foreach my $srcKey (keys %{$stateGraph}) {
    my $src = $self->cState ($srcKey);

    # only Null states
    unless ($self->isNull ($src)) { next; }

    my $statesDone = {};
    $nullConnections->{$srcKey} = [];
    $self->buildNullConnectionsHelper ($src, $nullConnections->{$srcKey}, $statesDone, $stateGraph);
  }

  return $nullConnections;
}

# Recursive helper for buildNullConnections().
# statesDone is used to prevent an infinite recursion in case of
# null cycles.
sub buildNullConnectionsHelper {
  my ($self, $beginState, $connectedStates, $statesDone, $stateGraph) = @_;

  # sanity check: only null states
  unless ($self->isNull ($beginState)) { croak ($self->cStateKey ($beginState), "is not a null state!\n"); }

  # if we've been here before, then return (prevent infinite recursion)
  if (defined $statesDone->{$beginState}) { return; }

  # update bookkeeping
  $statesDone->{$beginState} = 1;

  foreach my $destKey (@{$stateGraph->{$self->cStateKey ($beginState)}}) {
    my $dest = $self->cState ($destKey);

    # store non-Null states
    if (!$self->isNull ($dest)) { push (@$connectedStates, $dest); next; }
    
    # call recursively
    $self->buildNullConnectionsHelper ($dest, $connectedStates, $statesDone, $stateGraph);
  }

}

# Returns an array of state keys of Null states which can transition directly to End.
sub buildNullToEnd {
  my ($self, $stateGraph) = @_;
  
  my $states = [];
  foreach my $srcKey (keys %{$stateGraph}) {
    my $src = $self->cState ($srcKey);

    # only Null states
    unless ($self->isNull ($src)) { next; }

    foreach my $destKey (@{$stateGraph->{$src}}) {
      my $dest = $self->cState ($destKey);
      
      # store current if connected to End
      if ($self->isEnd ($dest)) { push (@$states, $srcKey); last; }
    }

  }

  return $states;
}

# Recursive helper for buildNullToEndConnections().
# statesDone is used to prevent an infinite recursion in case of
# null cycles.
sub buildNullToEndHelper {
  my ($self, $beginState, $states, $statesDone, $stateGraph) = @_;

  # sanity check: only null states
  unless ($self->isNull ($beginState)) { croak ($self->cStateKey ($beginState), "is not a null state!\n"); }

  # if we've been here before, then return (prevent infinite recursion)
  if (defined $statesDone->{$beginState}) { return; }

  # update bookkeeping
  $statesDone->{$beginState} = 1;

  foreach my $destKey (@{$stateGraph->{$self->cStateKey ($beginState)}}) {
    my $dest = $self->cState ($destKey);

    # store current if connected to End
    if ($self->isEnd ($dest)) { push (@$states, $beginState); next; }
    
    # call recursively
    $self->buildNullToEnd ($dest, $states, $statesDone, $stateGraph);
  }

}


# Constructs the reduced state graph (Null states removed, 
# excepting those which are children of bifurcation states).
sub buildReducedStateGraph {
  my ($self) = @_;

  my $reducedStateGraph = {};
  my $nullConnections = $self->buildNullConnections ($self->stateGraph);

  # construct list of Null states which are children of bifurcation states
  my $okNullStates = {};
  foreach my $srcKey (keys %{$self->stateGraph}) {
    my $src = $self->cState ($srcKey);
    unless ($self->isBifurc ($src)) { next; }
    my $children = $self->bifurcComposedChildren ($src);
    foreach my $branch (keys %{$children}) {
      my $child = $children->{$branch};
      if ($self->isNull ($child)) { $okNullStates->{$self->cStateKey ($child)} = 1; }
    }
  }

  foreach my $srcKey (keys %{$self->stateGraph}) {
    # skip null states (meaning don't add them to the reducedStateGraph) 
    # unless they're children of bifurcation states
    if ($self->isNull ($self->cState ($srcKey))) {
      unless ($okNullStates->{$srcKey}) { next; }
    }

    # leave bifurcation states as in the original stateGraph
    if ($self->isBifurc ($self->cState ($srcKey))) {
      foreach my $destKey (@{$self->stateGraph->{$srcKey}}) {
	push (@{$reducedStateGraph->{$srcKey}}, $destKey);
      }
      next;
    }

    # store dummy array reference for End state (just to be consistent with stateGraph)
    if ($self->isEnd ($self->cState ($srcKey))) {
      $reducedStateGraph->{$srcKey} = [qw //];
    }

    foreach my $destKey (@{$self->stateGraph->{$srcKey}}) {
      # if dest is a null state, then add transitions as appropriate
      if ($self->isNull ($self->cState ($destKey))) {
	# add a transition to all non-Null states connected by Null states
	if (ref $nullConnections->{$destKey}) {
	  foreach my $newDest (@{$nullConnections->{$destKey}}) {
	    my $newDestKey = $self->cStateKey ($newDest);
	    push (@{$reducedStateGraph->{$srcKey}}, $newDestKey);
	  }
	} else { croak "Oops, no nullConnections entry for '$destKey'\n"; }
      }
      # if dest isn't a null state, then store the transition
      else {
	push (@{$reducedStateGraph->{$srcKey}}, $destKey);
      }

    }
  }

  return $reducedStateGraph;
}


sub buildFullyReducedBifurcMap {
  my ($self) = @_;

  my $nullConnections = $self->buildNullConnections ($self->reducedStateGraph);

  my $reducedMap = {};
  foreach my $srcKey (keys %{$self->bifurcMap}) {
    # get the original children of the bifurcation
    my $children = $self->bifurcComposedChildren ($self->cState ($srcKey));

    # loop over branches (left, center, right)
    foreach my $branch (keys %{$children}) {
      my $child = $children->{$branch};
      my $childKey = $self->cStateKey ($child);

      # if the child is Null, then add all accessible states 
      # as bifurcation child states
      if ($self->isNull ($child)) {
	if (ref $nullConnections->{$childKey}) {
	  foreach my $newChild (@{$nullConnections->{$childKey}}) {
	    my $newChildKey = $self->cStateKey ($newChild);
	    $reducedMap->{$srcKey}->{$branch}->{$newChildKey} = "(" . $self->effTransProb ($srcKey, $childKey) . ") * (" . $self->effTransProb ($childKey, $newChildKey) . ")";
	    if ($self->directives_debug) {
	      print "Added $srcKey -> $newChildKey = ", $reducedMap->{$srcKey}->{$branch}->{$newChildKey}, "\n";
	    }
	  }
	} else { croak "Error, no nullConnections entry for '$childKey'\n"; }
      } else {
	$reducedMap->{$srcKey}->{$branch}->{$childKey} = $self->bifurcMap->{$srcKey}->{$branch}->{$childKey};
      }
    }

  }

  return $reducedMap;
}

# Constructs the fully reduced state graph (Null states removed).
# Built using the reducedStateGraph.
sub buildFullyReducedStateGraph {
  my ($self) = @_;

  my $fullyReducedStateGraph = {};
  my $nullConnections = $self->buildNullConnections ($self->reducedStateGraph);

  foreach my $srcKey (keys %{$self->reducedStateGraph}) {
    # skip null states
    if ($self->isNull ($self->cState ($srcKey))) { next; }

    # store dummy array reference for End state (just to be consistent with stateGraph)
    if ($self->isEnd ($self->cState ($srcKey))) {
      $fullyReducedStateGraph->{$srcKey} = [qw //];
    }

    foreach my $destKey (@{$self->reducedStateGraph->{$srcKey}}) {
      # if dest is a null state, then add transitions as appropriate
      if ($self->isNull ($self->cState ($destKey))) {
	# add a transition to all non-Null states connected by Null states
	if (ref $nullConnections->{$destKey}) {
	  foreach my $newDest (@{$nullConnections->{$destKey}}) {
	    my $newDestKey = $self->cStateKey ($newDest);
	    push (@{$fullyReducedStateGraph->{$srcKey}}, $newDestKey);
	  }
	} else { croak "Oops, no nullConnections entry for '$destKey'\n"; }
      }
      # if dest isn't a null state, then store the transition
      else {
	push (@{$fullyReducedStateGraph->{$srcKey}}, $destKey);
      }

    }
  }

  return $fullyReducedStateGraph;
}



################
# Show methods #
################


sub show {
  my ($self) = @_;

  # state graph and tm
  if ($self->directives_graph || $self->directives_all) {
    $self->showStateGraph ("stateGraph", $self->stateGraph);
  }
  if ($self->directives_tm || $self->directives_all) {
    print $self->comment, " Transition matrix\n";
    $self->showTM ($self->tm);
    $self->showBifurcMap ($self->bifurcMap);
  }

  if ($self->directives_rgraph || $self->directives_all) {
    $self->showStateGraph ("reducedStateGraph", $self->reducedStateGraph);
  }
  if ($self->directives_rtm || $self->directives_all) {
    print $self->comment, " Reduced transition matrix\n";
    $self->showTM ($self->reducedTM);
    $self->showBifurcMap ($self->bifurcMap);
  }

  if ($self->directives_frgraph || $self->directives_all) {
    $self->showStateGraph ("fullyReducedStateGraph", $self->fullyReducedStateGraph);
  }
  if ($self->directives_frtm || $self->directives_all) {
    print $self->comment, " Fully-reduced transition matrix\n";
    $self->showTM ($self->fullyReducedTM);
    $self->showBifurcMap ($self->fullyReducedBifurcMap);
  }

  # don't show redundant states or statistics if creating GraphViz output
  unless ($self->directives_graph || $self->directives_rgraph || $self->directives_frgraph) {
    # redundant states
    if ($self->directives_redundant) { $self->showRedundant ($self->tm); }

    # statistics
    print "\n", $self->comment, " Singlet machine:\n";
    $self->singlet->showStatistics();
    print $self->comment, " Branch machine:\n";
    $self->branch->showStatistics();
    print $self->comment, " Composed machine:\n";
    $self->showStatistics ($self->tm, $self->bifurcMap);
    print $self->comment, " Reduced composed machine:\n";
    $self->showStatistics ($self->reducedTM, $self->bifurcMap);
    print $self->comment, " Fully-reduced composed machine:\n";
    $self->showStatistics ($self->fullyReducedTM, $self->fullyReducedBifurcMap);
  }

}

1
