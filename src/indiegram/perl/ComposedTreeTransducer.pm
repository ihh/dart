#!/usr/bin/perl -w

package ComposedTreeTransducer;

use Carp;
use strict;
use Util;

use vars '@ISA';

sub new {
  my ($class, $directives, $comment) = @_;

  my $self = {
	      'directives' => {},
	      'startState' => {}, # the "overall" start state of the model
	      'endState' => {},
	      'stateGraph' => {},
	      'tm' => {},
	      'bifurcMap' => {},
	      'comment' => "" # the "comment" delimiter (e.g. '//' for C++) (defaults to ';')
	     };

  bless $self, ref ($class) || $class;

  $self->directives ($directives);
  if (defined $comment) { $self->comment ($comment); }
  else { $self->comment (";"); }

  $self->_initialize();

  return $self;
}


# catch methods
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD; # $AUTOLOAD contains the fully qualified name of the original subroutine
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq 'DESTROY';

  # check for directives accessor, e.g. $self->directives_debug or $self->directives_('debug')
  if ($sub =~ /^directives_(\S*)$/i) {
    my $flag = lc ($1);
    $flag = shift @args unless (length $flag); # this catches the argument 'debug' in the second example usage given above
    if (!defined $self->{'directives'}->{$flag}) {
      $self->{'directives'}->{$flag} = "";     # if no such flag exists, create one with the value ""
    }                                          # we therefore have to test 'if ($self->directives_debug)' rather than 'if (defined $self->directives_debug)'
    return $self->{'directives'}->{$flag};     # the second will always be true because AUTOLOAD will create an empty flag for us
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

sub _initialize {
  my ($self) = @_;

}

##########
# Basics #
##########

# A composed state, %cState, is a hash mapping a node to
# the state of the TT (or SCFG) living on its parent branch.
# The corresponding state key, $cStateKey,
# is an underscore-delimited string of states.
# Given a %cState, return the corresponding $cStateKey.
sub cStateKey {
    my ($self, $cState) = @_;
    # map BLOCK LIST
    #  Evaluates EXPR or BLOCK for each element of the LIST, locally setting $_ to refer to the element. Modifying $_ will modify the corresponding element from LIST. Returns the list of results.
#    my @branchStates = map ($cState->{$_}, 1..$self->nodes);
    my @branchStates = map ($cState->{$_}, 1..scalar(keys %$cState));
    return join ("_", map (defined ($_) ? $_ : "", @branchStates));
}

# Given a $cStateKey, return the corresponding %cState.
sub cState {
    my ($self, $cStateKey) = @_;
    my @branchStates = split (/_/, $cStateKey);
    my $cState;
    for (my $i = 0; $i < scalar (@branchStates); ++$i) {
	if ($branchStates[$i] ne "") {
	    $cState->{$i+1} = $branchStates[$i]; # +1 indexing so that properly keyed by node
	}
    }
    return $cState;
}


######################################
# Methods which must be implemented! #
######################################

sub stateCmp;

sub isStart;
sub isEnd;
sub isEmit;
sub isNull;
sub isBifurc;

sub areEquivalent;

########################
# State access methods #
########################


####################
# Graph algorithms #
####################

# Find redundant states (states with identical rows in the TM, emit
# profiles and distributions) in the passed state graph.
# Can be called on e.g. any of the defined state graphs in FourWayComposedTT,
# including e.g. stateGraph, reducedStateGraph and fullyReducedStateGraph.
sub findRedundant {
  my ($self, $tm) = @_;

  my $redundant = {};
  
  foreach my $srcKey (keys %$tm) {
    foreach my $candKey (keys %$tm) {
      unless ($srcKey ne $candKey) { next; }
      if ($self->areEquivalent ($srcKey, $candKey, $tm)) {
	unless (exists $redundant->{$candKey} and
		defined $redundant->{$candKey}->{$srcKey})
	  { $redundant->{$srcKey}->{$candKey} = 1; }
      }
      else { next; }
    }
  }

  return $redundant;
}

# to do: convert all this to use a TM and bifurcMap 
# (b/c this is how I'll have to do it for Indiegram)

# Given a state graph, perform a topological sort.
# Warns if a cycle found, but otherwise is NOT robust w.r.t. cycles:
# doesn't sort nodes encountered after finding a cycle.
sub topologicalSort {
  my ($self, $stateGraph) = @_;

  # initialize and fill sources and sinks
  my $sources = {}; # $sources->{$sKey} is the number of parents of $sKey
  my $sinks = {}; # $sinks->{$sKey} is an array of child states of $sKey
  foreach my $srcKey (keys %{$stateGraph}) {
    $sources->{$srcKey} = 0;
    $sinks->{$srcKey} = [];
  }
  foreach my $srcKey (keys %{$stateGraph}) {
    foreach my $destKey (@{$stateGraph->{$srcKey}}) {
      $sources->{$destKey} += 1;
      push (@{$sinks->{$srcKey}}, $destKey);
    }
  }

  # initialize unsorted and sorted arrays
  my @unsorted;
  my @sorted;
  foreach my $sKey (keys %{$stateGraph}) {
    push (@unsorted, $sKey);
  }

  # do a topological sort
  while (scalar (@unsorted)) {

    # look for nodes with no parents
    my $found = 0;
    my $node;
    foreach my $sKey (@unsorted) {
      if ($sources->{$sKey} == 0) { # if $sKey has no parents
	$node = $sKey;
	$found = 1;
	last;
      }
    }

    # if no parentless node found, then exit loop
    if (!$found) { last; }

    # store the node in the sorted array
    push (@sorted, $node);

    # remove the stored node from the graph (update sources and sinks)
    Util->aRemove (\@unsorted, $node);
    foreach my $destKey (@{$sinks->{$node}}) {
      $sources->{$destKey} -= 1;
    }
  }

  if (scalar (@unsorted)) {
    print $self->comment, " WARNING: Topological sort failed due to cycles.  ", scalar (@unsorted), " unsorted nodes (of ", scalar (@unsorted) + scalar (@sorted)," nodes total).\n";
    push (@sorted, @unsorted);
  }
  
  return \@sorted;
}

# Given a state graph, perform a pseudo-topological sort.
# More strictly, collapses each strongly connected component into a single node
# and then topologically sorts the resulting graph.
sub pseudoTopologicalSort {
  my ($self, $stateGraph) = @_;

  # find all component of the graph
  # collapse each into a single node
  # topologically sort the resulting DAG
  # sort each component
  # insert the sorted component nodes into the correct
  # location in the topological sort

}

# Given a state graph and a strongly connected component of the graph,
# sort the nodes in the component by backtracking from the 
# "exit" states of the component.
sub sortComponent {
  my ($self, $component, $stateGraph) = @_;

  

}

# Given a state graph, identify all strongly connected components.
#  (strongly connected components = maximal strongly connected subgraphs)
sub findConnectedComponents {
  my ($self, $stateGraph) = @_;

  # build list of connected components
  my $dfsnum; # this must be a /reference/ to a scalar b/c the dfsnum of the DFS must be globally visible
  $$dfsnum = 0;
  my $dfsnumMap = {}; # the DFS index for each node in the search tree
  my $lowlinkMap = {};
  my $stack = [];
  my $components = [];
  foreach my $sKey (keys %{$stateGraph}) {
    unless (defined $dfsnumMap->{$sKey}) {
      $self->tarjanHelper ($sKey, $dfsnum, $dfsnumMap, $lowlinkMap, $stack, $components, $stateGraph);
    }
  }

  return $components;
}

# Recursive helper for findConnectedComponents.
# This is Tarjan's algorithm:
# see http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm.
#Input: Graph G = (V, E), Start node v0
#
#index = 0                       // DFS node number counter 
#S = empty                       // An empty stack of nodes
#tarjan(v0)                      // Start a DFS at the start node
#
#procedure tarjan(v)
#  v.index = index               // Set the depth index for v
#  v.lowlink = index
#  index = index + 1
#  S.push(v)                     // Push v on the stack
#  forall (v, v') in E do        // Consider successors of v 
#    if (v'.index is undefined)  // Was successor v' visited? 
#      tarjan(v')                // Recurse
#      v.lowlink = min(v.lowlink, v'.lowlink)
#    elseif (v' in S)            // Is v' on the stack?
#      v.lowlink = min(v.lowlink, v'.index)
#  if (v.lowlink == v.index)     // Is v the root of an SCC?
#    print "SCC:"
#    repeat
#      v' = S.pop
#      print v'
#    until (v' == v)
sub tarjanHelper {
  my ($self, $sKey, $dfsnum, $dfsnumMap, $lowlinkMap, $stack, $components, $stateGraph) = @_;

  # push the current node
  push (@$stack, $sKey);

  # store and increment indices
  $dfsnumMap->{$sKey} = $$dfsnum;
  $lowlinkMap->{$sKey} = $$dfsnum;
  $$dfsnum += 1;

  # DFS
  foreach my $dKey (@{$stateGraph->{$sKey}}) {
    # if sKey->dKey is a tree arc (== dKey not yet visited in the DFS)
    if (!defined $dfsnumMap->{$dKey}) {
      $self->tarjanHelper ($dKey, $dfsnum, $dfsnumMap, $lowlinkMap, $stack, $components, $stateGraph);
      $lowlinkMap->{$sKey} = Util->min ($lowlinkMap->{$sKey}, $lowlinkMap->{$dKey});
    }
    # if sKey->dKey is a back-edge or cross-link
    # NB: The condition that dKey be on the stack is equivalent to requiring that
    # (dfsnum (dKey) < dfsnum (sKey)) && (dKey in the same component as sKey)
    #    (if it were in a different component then it would already have been popped from the stack)
    elsif (Util->aContains ($stack, $dKey)) {
      $lowlinkMap->{$sKey} = Util->min ($lowlinkMap->{$sKey}, $dfsnumMap->{$dKey});
    }
  }

  # if at the root of a subtree, then found a strongly connected component
  if ($lowlinkMap->{$sKey} == $dfsnumMap->{$sKey}) {
    my @path;
    my $dKey = "";
    while ($dKey ne $sKey) {
      $dKey = pop (@$stack);
      push (@path, $dKey);
    }
    push (@$components, \@path);
  }

}



################
# Show methods #
################


# Display the state graph of the composed machine
# in GraphViz format.
# Sanitizes the names of bifurc states for GraphViz output.
sub showStateGraph {
  my ($self, $name, $stateGraph) = @_;

  print "digraph $name \{\n";
  # transitions from non-Bifurc sorted states
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$stateGraph}) { # $srcKey is a composed state key
    if ($self->isBifurc ($self->cState ($srcKey))) { next; }
    foreach my $destKey (sort {$self->stateCmp ($a, $b)} @{$stateGraph->{$srcKey}}) {
      print "\t", Util->toAlphanumeric ($srcKey), " -> ", Util->toAlphanumeric ($destKey), ";\n";
    }
    print "\n";
  }

  # "transitions" from bifurcations
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$stateGraph}) { # $srcKey is a composed state key
    unless ($self->isBifurc ($self->cState ($srcKey))) { next; }
    # special attributes for bifurcation nodes
    print Util->toAlphanumeric ($srcKey), " \[shape=box, style=filled\];\n";
    foreach my $destKey (sort {$self->stateCmp ($a, $b)} @{$stateGraph->{$srcKey}}) {
      if (!$self->directives_nosuppress && $self->isEnd ($self->cState ($destKey))) { next; }
      print "\t", Util->toAlphanumeric ($srcKey), " -> ", Util->toAlphanumeric ($destKey), " \[style=bold\];\n";
    }
    print "\n";
  }

  print "\}\n";
}

# Display the transition matrix of the composed machine.
sub showTM {
  my ($self, $tm) = @_;

  # display sorted states
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$tm}) {
    my $src = $self->cState ($srcKey);
      
    print "; $srcKey (";
    if ($self->isStart ($src)) {
      print "Start";
    } elsif ($self->isEnd ($src)) {
      print "End";
    } elsif ($self->isEmit ($src)) {
      print "Emit: ", join (",", $self->emitSet ($src));
    } elsif ($self->isNull ($src)) { 
      print "Null";
    } elsif ($self->isBifurc ($src)) {
      croak "Uh-oh!  There shouldn't be bifurcs in the TM, but I encountered '$srcKey'.\n";
    } else {
      croak "Unreachable code.\n";
    }
    print ")\n";
    foreach my $destKey (sort {$self->stateCmp ($a, $b)} keys %{$tm->{$srcKey}}) {
      print "   $srcKey -> $destKey  {", $tm->{$srcKey}->{$destKey}, "};\n";
    }
  }

}

sub showBifurcMap {
  my ($self, $bifurcMap) = @_;

  if (scalar keys %{$bifurcMap}) {
    print "\n;; Bifurcations\n";
  } else { print "\n;; No bifurcations.\n"; }
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$bifurcMap}) {
    print "; $srcKey (Bifurc)\n";
    # loop over left, center, right children
    foreach my $branch (keys %{$bifurcMap->{$srcKey}}) {
      # don't display branches with only End transitions unless specifically requested
      my $cntEnd = 0;
      foreach my $childKey (keys %{$bifurcMap->{$srcKey}->{$branch}}) {
	if ($self->isEnd ($self->cState ($childKey))) { $cntEnd += 1; }
      }
      unless ($self->directives_nosuppress) {
	if ($cntEnd eq scalar (keys %{$bifurcMap->{$srcKey}->{$branch}})) { next; }
      }
      print ";  ($branch):\n";
      foreach my $childKey (keys %{$bifurcMap->{$srcKey}->{$branch}}) {
	print "      $srcKey -> $childKey \{", $bifurcMap->{$srcKey}->{$branch}->{$childKey}, "\};\n";
      }
    }
    print "\n";
  }

}

sub showRedundant {
  my ($self, $tm) = @_;

  my $redundant = $self->findRedundant ($tm);
  if (scalar keys %{$redundant}) {
    print $self->comment, " Redundant states in the reduced TM:\n";
  } else { print "\n", $self->comment, " No redundant states in the reduced TM.\n"; }
  while (my ($s1, $s2hash) = each %{$redundant}) {
    while (my ($s2, $discard) = each %{$s2hash}) {
      print $self->comment, "  $s1 <=> $s2\n";
    }
  }

}

# Show statistics on the composed machine.
sub showStatistics {
  my ($self, $tm, $bifurcMap) = @_;

  my $numTMStates = scalar (keys %{$tm});

  # count the number of Null and Emit states;
  my $numEmit = 0;
  my $numNull = 0;
  foreach my $sKey (keys %{$tm}) {
    my $s = $self->cState ($sKey);
    if ($self->isEmit ($s)) { $numEmit += 1; }
    if ($self->isNull ($s)) { $numNull += 1; }
  }

  # count the number of transitions in the TM
  my $numTrans = 0;
  foreach my $srcKey (keys %{$tm}) {
    foreach my $destKey (keys %{$tm->{$srcKey}}) {
      $numTrans += 1;
    }
  }
  # count the number of non-deterministic bifurcations
  my $numDeterministicBifurc = 0;
  foreach my $srcKey (keys %{$bifurcMap}) {
    my $tmp = 1;
    foreach my $branch (keys %{$bifurcMap->{$srcKey}}) {
      $tmp *= scalar (keys %{$bifurcMap->{$srcKey}->{$branch}});
    }
    $numDeterministicBifurc += $tmp;
  }

  # NB: +1 to account for the End state (Start state is already present in the TM as a source state)
  my $numStates = $numTMStates + $numDeterministicBifurc + 1;
  print $self->comment, "   $numStates = $numTMStates+$numDeterministicBifurc+1 effective total states\n";
  print $self->comment, "   $numDeterministicBifurc Bifurcation states (effectively deterministic)\n";
  print $self->comment, "      (", scalar (keys %{$bifurcMap}), " possibly non-deterministic bifurc states)\n";
  print $self->comment, "   $numEmit Emit states\n";
  print $self->comment, "   $numNull Null states\n";
  print $self->comment, "   $numTMStates (non-bifurcation) states in TM\n";
  print $self->comment, "   $numTrans transitions in TM\n";
  print $self->comment, "      $numTMStates^2 = ", scalar (keys %{$tm})*scalar (keys %{$tm}), " (for comparison)\n";
}


1
