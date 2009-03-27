#!/usr/bin/perl -w

package ComposedTreeTransducer::TripletSCFG;

use strict;
use Carp;
use Util;

use ComposedTreeTransducer;
use ComposedTreeTransducer::FourWayComposedTT;

our @ISA = qw /ComposedTreeTransducer/;
use vars '@ISA';

# To do:
#  - 

# This is designed specifically for use with my indiegram/tripletscfg code.
# There's no fiddly stuff with transducer composition or the like.  State
# typing is simple: Start, Null, Emit, Bifurc and End.

my @seqs = qw /X Y Z/;

sub new {
  my ($class, $singletfile, $branchfile, $directives) = @_;

  my $self = ComposedTreeTransducer->new ($directives, "//");

  bless $self, ref ($class) || $class;

  # store the "parent" FourWayComposedTT
  $self->{'fourWay'} = "";
  $self->fourWay (ComposedTreeTransducer::FourWayComposedTT->new ($singletfile, $branchfile));
  
  $self->_initialize();

  return $self;
}

sub _initialize {
  my ($self) = @_;
  
  # store (reduced) state graph
  $self->stateGraph ($self->fourWay->reducedStateGraph);

  # store (reduced) TM
  # parse trans prob params into the format we want
  foreach my $srcKey (keys %{$self->fourWay->reducedTM}) {
    foreach my $destKey (keys %{$self->fourWay->reducedTM->{$srcKey}}) {
      $self->tm->{$srcKey}->{$destKey} = $self->formatParams ($self->fourWay->reducedTM->{$srcKey}->{$destKey});
    }
  }
  # hack in an "overall" Start state
  # NB: This is distinct from the startState of a FourWayComposedTT.  That startState is
  # a composed state composed of the start states of the singlet and branch machines and can be
  # a child of a bifurcation; e.g., it's not really distinct from the other states of type Start in the model.
  # This startState, in contrast, can be reached from no other state of the model.
  # The Triplet_SCFG code assumes that the model has such an overall startState (-1 per Grammar_state_enum).
  # Deterministic transition to the startState of the corresponding FourWayComposedTT.
  $self->startState ($self->cState ('overallStart_overallStart_overallStart_overallStart'));
  $self->tm->{$self->cStateKey ($self->startState)}->{$self->cStateKey ($self->fourWay->startState)} = '(1)';

  # store bifurcMap
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$self->fourWay->bifurcMap}) {
    foreach my $branch (keys %{$self->fourWay->bifurcMap->{$srcKey}}) {
      # sanity check: bifurcations had better be deterministic in the reduced bifurcMap!
      if (scalar (keys %{$self->fourWay->bifurcMap->{$srcKey}->{$branch}}) ne 1) {
	croak ("$branch branch of bifurcation state '$srcKey' must have exactly 1 child.");
      }
      my ($childKey) = keys %{$self->fourWay->bifurcMap->{$srcKey}->{$branch}};
      # only left-bifurcations allowed
      if ($branch eq 'right') {
	if (!$self->isEnd ($self->cState ($childKey))) {
	  croak ("only left-bifurcations allowed: bifurcation state '$srcKey' has a non-End right child '$childKey'.");
	}
	next;
      }
      $self->bifurcMap->{$srcKey}->{$branch}->{$childKey} = $self->fourWay->bifurcMap->{$srcKey}->{$branch}->{$childKey};
    }
  }

  # find all strongly connected components of the graph on the Null and Bifurc states
  # create lists of all Null and Bifurc states
  my $nullStates = {};
  foreach my $sKey (keys %{$self->tm}) {
    if ($self->isNull ($self->cState ($sKey))) {
      $nullStates->{$sKey} = 1;
    }
  }
  my $bifurcStates = {};
  while (my ($sKey, $discard) = each %{$self->bifurcMap}) {
    $bifurcStates->{$sKey} = 1;
  }
  # construct a state graph which contains only the desired states
  my $tmpGraph = {};
  foreach my $sKey (keys %{$self->stateGraph}) {
    if (defined $nullStates->{$sKey} || defined $bifurcStates->{$sKey}) {
      foreach my $dKey (@{$self->stateGraph->{$sKey}}) {
	if (defined $nullStates->{$dKey} || defined $bifurcStates->{$dKey}) {
	  push (@{$tmpGraph->{$sKey}}, $dKey);
	}
      }
    }
  }

  # get the components
  my $components = $self->findConnectedComponents ($tmpGraph);

  # to do: eventually use them when want to do Inside/Outside instead of just CYK

}

# Get the root sequence name.
sub root {
  my ($class) = @_;
  return 'W';
}

# Map nodes to sequence names.
sub seq {
  my ($class, $n) = @_;
  my %map = (
	     1 => 'W',
	     2 => 'X',
	     3 => 'Y',
	     4 => 'Z'
	    );
  if (exists $map{$n}) { return $map{$n}; }
  else { croak "No entry for node '$n'."; }
}

# Return a hash giving the emit profile for W, X, Y, Z.
sub emitProfile {
  my ($self, $s) = @_;

  my $emitProfile = {};
  my $tmp = $self->fourWay->emitProfile ($s);
  while (my ($n, $p) = each %{$tmp}) {
    $emitProfile->{$self->seq ($n)} = $p;
  }

  return $emitProfile;
}

# Return a hash giving the emit prob factor for W, X, Y, Z.
# Replaces e.g. m_2(lr_2|lr_1) with m_x(xlr|wlr).
sub emitProbs {
  my ($self, $s) = @_;

  my $emitProbs = {};
  my $tmp = $self->fourWay->emitProbs ($s);
  while (my ($n, $p) = each %{$tmp}) {
    # replace e.g. m_2(lr_2|lr_1) with m_2(lr_1,lr_2)
    $p =~ s/\((\w+)\|(\w+)\)/\($2,$1\)/;
    # replace e.g. lr_2 with xlr
    $p =~ s/(lr)_([1-4])/lc ($self->seq ($2)) . 'lr'/ge;
    $p =~ s/([l,r])_([1-4])/lc ($self->seq ($2)) . $1/ge;
    # replace e.g. m_2(xl|wl) with m_x(xl|wl)
    $p =~ s/_([1-4])/'_' . lc ($self->seq ($1))/e;
    $emitProbs->{$self->seq ($n)} = $p;
  }

  return $emitProbs;
}

# Is the unique overall Start state specified by $self->startState?
sub isStart {
  my ($self, $s) = @_;
  unless ($self->cStateKey ($self->startState) eq $self->cStateKey ($s)) { return 0; }
  return 1;
}

# Is the unique overall End state?
sub isEnd {
  my ($self, $s) = @_;
  if ($self->isStart ($s)) { return 0; } # handle this separately b/c not a state of fourWay
  return $self->fourWay->isEnd ($s);
}

# Is an Emit state?
# Only if there are emissions to an extant sequence.
sub isEmit {
  my ($self, $s) = @_;
  if ($self->isStart ($s)) { return 0; } # handle this separately b/c not a state of fourWay
  my $emitProfile = $self->emitProfile ($s);
  my $is = 0;
  while (my ($seq, $p) = each %{$emitProfile}) {
    if ($seq eq $self->root()) { next; }
    if ($p ne "") { $is = 1; }
  }
  return $is;
}

# Is it a Bifurc state?
# This is the same as for FourWayComposedTT.
sub isBifurc {
  my ($self, $s) = @_;
  if ($self->isStart ($s)) { return 0; } # handle this separately b/c not a state of fourWay
  return $self->fourWay->isBifurc ($s);
}

# Is it a Null state?
sub isNull {
  my ($self, $s) = @_;
  return (!($self->isStart ($s) || $self->isEnd ($s) || $self->isEmit ($s) || $self->isBifurc ($s)));
}

sub stateType {
  my ($self, $sKey) = @_;
  my $s = $self->cState ($sKey);

  my $t = "";
  if ($self->isStart ($s)) { $t = 'Start'; }
  elsif ($self->isEnd ($s)) { $t = 'End'; }
  elsif ($self->isNull ($s)) { $t = 'Null'; }
  elsif ($self->isBifurc ($s)) { $t = 'Bifurc'; }
  elsif ($self->isEmit ($s)) {
    $t = 'Emit';
    my $emitProfile = $self->emitProfile ($s);
    foreach my $seq (sort keys %{$emitProfile}) {
      if ($seq eq $self->root()) { next; }
      if ($emitProfile->{$seq} ne "") {
	$t .= uc ($seq . $emitProfile->{$seq});
      }
    }
    if ($t eq 'Emit') { croak "Unreachable code: state '$sKey'!\n"; }
  }
  else { croak "Unreachable code: state '$sKey'.\n"; }

  return $t;
}


# Sort by state type, then if both Emit sort alphabetically by state type.
# Break ties by sorting alphabetically by state key.
sub stateCmp {
  my ($self, $aKey, $bKey) = @_;

  my %weights = (
		 'Start' => 1,
		 'Null' => 10,
		 'Emit' => 100,
		 'Bifurc' => 1000,
		 'End' => 10000,
		 );
  
  my $aType = $self->stateType ($aKey);
  my $bType = $self->stateType ($bKey);
  # if not type Emit
  if (defined $weights{$aType} && defined $weights{$bType}) {
    if ($aType ne $bType) { return ($weights{$aType} <=> $weights{$bType}); }
    else { return $aKey cmp $bKey; } }
  elsif (defined $weights{$aType} && !defined $weights{$bType}) {
    return ($weights{$aType} <=> $weights{'Emit'}); }
  elsif (!defined $weights{$aType} && defined $weights{$bType}) {
    return ($weights{'Emit'} <=> $weights{$bType}); }
  else {
    if ($aType ne $bType) { return ($aType cmp $bType); }
    else { return ($aKey cmp $bKey); } }

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

  # now check emission probs for Emit states
  if ($t1 =~ /Emit/) {
    my $emitProb1 = $self->emitProbs ($cState1);
    my $emitProb2 = $self->emitProbs ($cState2);
    while (my ($seq, $p) = each %$emitProb1) {
      unless (exists $emitProb2->{$seq} and
	      $p eq $emitProb2->{$seq})
	{ return 0; }
    }
    while (my ($seq, $p) = each %$emitProb2) {
      unless (exists $emitProb1->{$seq} and
	      $p eq $emitProb1->{$seq})
	{ return 0; }
    }
  }

  return 1;
}

# Format the parameters to accomodate my Triplet_SCFG code.
# Currently does nothing.
sub formatParams {
  my ($self, $p) = @_;

  return $p;
}

# Given an array ref of state keys (all of same stateType) 
# write prettily-formatted state initialization and emissions stuff.
sub writeStateInit {
  my ($self, $states) = @_;
  
  my $t = $self->stateType ($states->[0]);

  print $self->comment, " state initialization: $t\n",
    $self->comment, "  (Automatically generated by 'tripletSCFG.pl' --write.)\n";
  foreach my $sKey (@$states) {
    print "init_emit ($sKey, $t);\n";
  }

}

# Given an array ref of state keys (all of same stateType),
# write prettily-formatted state emissions stuff.
# Currently can ONLY handle unpaired (left) and paired emissions.
sub writeStateEmit {
  my ($self, $states) = @_;
  
  my $t = $self->stateType ($states->[0]);
  my %c;

  my $isPaired = 0;

  # keep track of which indices we need to initialize over
  # unpaired (left) emissions
  if ($t =~ /^Emit(XL)?(YL)?(ZL)?$/) {
    if (defined $1) { $c{xl} = 1; }
    if (defined $2) { $c{yl} = 1; }
    if (defined $3) { $c{zl} = 1; } }
  # paired emissions
  elsif ($t =~ /^Emit(XLR)?(YLR)?(ZLR)?$/) {
    $isPaired = 1;
    if (defined $1) { $c{xl} = $c{xr} = 1; }
    if (defined $2) { $c{yl} = $c{yr} = 1; }
    if (defined $3) { $c{zl} = $c{zr} = 1; } }
  else { croak "I don't know what to do with states of type $t\n"; }

  my $xl = (exists $c{xl}) ? 'xl' : '0';
  my $xr = (exists $c{xr}) ? 'xr' : '0';
  my $yl = (exists $c{yl}) ? 'yl' : '0';
  my $yr = (exists $c{yr}) ? 'yr' : '0';
  my $zl = (exists $c{zl}) ? 'zl' : '0';
  my $zr = (exists $c{zr}) ? 'zr' : '0';

  print $self->comment, " emissions: $t\n",
    $self->comment, "  (Automatically generated by 'tripletSCFG.pl' --write.)\n";

  # print for loops over xl, xr, etc. as appropriate
  foreach my $c (sort keys %c) {
    print "for (int $c = 0; $c < SCFG_alphabet_size; ++$c)\n";
  }
  print "  {
    int emit_idx = (*this).emit_idx ($t, $xl, $xr, $yl, $yr, $zl, $zr);
    Prob pr;
";
  if ($isPaired) {
    if (exists $c{xr}) { print "    int xlr = (*this).emit_idx (EmitXLR, xl, xr, 0, 0, 0, 0);\n"; }
    if (exists $c{yr}) { print "    int ylr = (*this).emit_idx (EmitYLR, 0, 0, yl, yr, 0, 0);\n"; }
    if (exists $c{zr}) { print "    int zlr = (*this).emit_idx (EmitZLR, 0, 0, 0, 0, zl, zr);\n"; }
  }

  # Now loop over all states in the passed array.
  foreach my $sKey (@$states) {

    my $s = $self->cState ($sKey);

    # get the emission probability factor
    my $factors = $self->emitProbs ($s);
    my @factors;
    foreach my $n (sort keys %{$factors}) {
      push @factors, $factors->{$n};
    }
    my $factor = join (" * ", @factors);

    # See whether marginalization over the ancestral node is necessary:
    # 2 possible cases for Emit states:
    # 1) Insertion at a leaf.
    # 2) Emission at the root and a match at one or more leaves.
    # Explicit marginilization (Felsenstein-esque) is necessary
    # in the 2nd case.
    my $emitCnt = scalar (keys %$factors); # number of emitting nodes

    my $emitProfile = $self->emitProfile ($s);
    my $profile = "";
    while (my ($seq, $p) = each %$emitProfile) {
      if ($p ne "") { $profile = $p; last; }
    }

    # Case 1:
    if (!exists $factors->{$self->root()} && $emitCnt == 1) {
      print "    emit\[$sKey\]\[emit_idx\] = Prob2Score ($factor);\n";
    }
    # Case 2:
    elsif (exists $factors->{$self->root()} && $emitCnt > 1) {
      print "    pr = 0;";
      my $wl = lc ($self->root()) . 'l';
      my $wr = lc ($self->root()) . 'r';
      if (!$isPaired) {
	print "
    for (int $wl = 0; $wl < SCFG_alphabet_size; ++$wl)
      pr += $factor;
"; }
      else {
	print "
    for (int $wl = 0; $wl < SCFG_alphabet_size; ++$wl)
      for (int $wr = 0; $wr < SCFG_alphabet_size; ++$wr)
        {
          // hack to create an emit_idx hash for (wl,wr)
          int wlr = (*this).emit_idx (EmitXLR, wl, wr, 0, 0, 0, 0);
          pr += $factor;
        }
"; }
      print "    emit\[$sKey\]\[emit_idx\] = Prob2Score (pr);\n";
    } else { croak "Unreachable code.\n"; }

  }
  print "  }\n";

}



sub writeFormatted {
  my ($self) = @_;

  # write statistics
  print $self->comment, " Triplet_SCFG statistics:\n";
  $self->showStatistics ($self->tm, $self->bifurcMap);
  print "\n";

  # write number of states
  # -1 is to not count the startState
  # note that this assumes that the End state is not present as a source
  # state in the TM!
  my $numStates = scalar (keys %{$self->tm}) - 1 + scalar (keys %{$self->bifurcMap});
  print "/// number of states (not counting Start and End)
static const int num_states = $numStates;\n";
  
  # write state enum
  print "
/// Triplet_SCFG state enumeration
/*
 * Automatically generated by 'tripletSCFG.pl' --write.
 */
";
  print "enum Triplet_states \{\n";
  my $cnt = 0;
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm}) {
    if ($self->isStart ($self->cState ($srcKey))) { next; }
    print "$srcKey = $cnt, ";
    if ($cnt % 4 == 3) { print "\n"; }
    $cnt += 1;
  }
  print "\n";
  foreach my $srcKey (keys %{$self->bifurcMap}) {
    $srcKey = Util->toAlphanumeric ($srcKey);
    print "$srcKey = $cnt, ";
    if ($cnt % 4 == 3) { print "\n"; }
    $cnt += 1;
  }
  print "\};\n\n";

  # write ancestral state information
  print "
/// Triplet_SCFG ancestral state information
/*
 * Automatically generated by 'tripletSCFG.pl' --write.
 */
";
  print "
enum State_type_ancestral \{ WNull = 0,
EmitWR = 1, EmitWL = 2, EmitWLR = 3,
WBifurc = 4 \};
\n\n";
  print "map<int, State_type_ancestral> state_type_ancestral_map;\n\n";
  print "/* Emit and Null states */\n";
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm}) {
    my $s = $self->cState ($srcKey);
    if ($self->isStart ($s)) { next; }
    my $prof = $self->emitProfile ($self->cState ($srcKey));
    if (exists $prof->{$self->root()}) {
      print "state_type_ancestral_map.insert (make_pair ($srcKey, SCFG_state_typing::Emit", uc ($self->root()), uc ($prof->{$self->root()}), "));\n";
    } elsif ($self->isBifurc ($s)) {
      die "Unreachable.\n";
    } else {
      print "state_type_ancestral_map.insert (make_pair ($srcKey, SCFG_state_typing::NullW));\n";
    }
  }
  print "\n/* Bifurc states */\n";
  foreach my $srcKey (keys %{$self->bifurcMap}) {
    my $s = $self->cState ($srcKey);
    if ($self->isStart ($s)) { next; }
    my $prof = $self->emitProfile ($self->cState ($srcKey));
    if ($self->isBifurc ($s)) {
      $srcKey = Util->toAlphanumeric ($srcKey);
      print "state_type_ancestral_map.insert (make_pair ($srcKey, SCFG_state_typing::BifurcW));\n";
    } else {
      die "Unreachable.\n";
    }
  }
  print "\n\n";

  # group Null and Emit states by stateType
  my @null;
  my $unpairedMap = {};
  my $pairedMap = {};
  foreach my $sKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm}) {
    my $t = $self->stateType ($sKey);
    if ($t eq 'Null') { push (@null, $sKey); }
    elsif ($t =~ 'Emit') {
      if ($t =~ /^Emit([X-Z]LR){1,3}$/) { push (@{$pairedMap->{$t}}, $sKey); }
      elsif ($t =~ /^Emit([X-Z][L,R]){1,3}$/) { push (@{$unpairedMap->{$t}}, $sKey); }
      else { croak "I don't know what to do with state '$sKey' (type $t)\n"; }
    }
  }

  # write state init: null
  print "/*\n * state initialization: Null states\n */\n";
  $self->writeStateInit (\@null);
  
  # write state init: unpaired
  print "\n/*\n * state initialization: unpaired Emit states\n */\n";
  foreach my $sKeys (sort keys %$unpairedMap) {
    $self->writeStateInit ($unpairedMap->{$sKeys});
    print "\n";
    $self->writeStateEmit ($unpairedMap->{$sKeys});
    print "\n";
  }

  # write state init: paired
  print "\n/*\n * state initialization: paired Emit states\n */\n";
  foreach my $sKeys (sort keys %$pairedMap) {
    $self->writeStateInit ($pairedMap->{$sKeys});
    print "\n";
    $self->writeStateEmit ($pairedMap->{$sKeys});
    print "\n";
  }

  
  # write state init: bifurcation
  print "/*\n * state initialization: Bifurc states\n */\n";
  foreach my $srcKey (keys %{$self->bifurcMap}) {
    my ($leftKey) = keys %{$self->bifurcMap->{$srcKey}->{left}};
    my ($centerKey) = keys %{$self->bifurcMap->{$srcKey}->{center}};
    $srcKey = Util->toAlphanumeric ($srcKey);
    $leftKey = Util->toAlphanumeric ($leftKey);
    $centerKey = Util->toAlphanumeric ($centerKey);
    print "init_bifurc ($srcKey, $leftKey, $centerKey);\n";
  }
  print "\n";

  # write transition_scores: all states!
  print "/*\n * transition_scores initialization\n */\n";
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm}) {
    my $s = $self->cState ($srcKey);
    if ($self->isStart ($s)) {
      foreach my $destKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm->{$srcKey}}) {
	my $d = $self->cState ($destKey);
	# check for End case
	if ($self->isEnd ($d)) { print "transition_scores.start_to_end() = "; }
	else { print "transition_scores.start\[", Util->toAlphanumeric ($destKey), "\] = "; }
	print "Prob2Score (", $self->tm->{$srcKey}->{$destKey}, ");\n";
      } }
    elsif ($self->isEnd ($s)) { next; }
    else {
      print "/* $srcKey */\n";
      foreach my $destKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm->{$srcKey}}) {
	my $d = $self->cState ($destKey);
	# sanity check: Start should never appear as a dest state
	if ($self->isStart ($d)) { croak ("'", $self->cStateKey ($self->startState), "' should never appear as a dest state!\n"); }
	# check for End case
	if ($self->isEnd ($d)) { print "transition_scores.end\[$srcKey\] = "; }
	else { print "transition_scores.transition ($srcKey, ", Util->toAlphanumeric ($destKey), ") = "; }
	print "Prob2Score (", $self->tm->{$srcKey}->{$destKey}, ");\n";
	} }
    print "\n";
  }

  # write transition_scores: effective direct transitions
  # ideally this would be handled directly within the FourWayComposedTransducer code,
  # but for technical reasons (because of how probabilities are stored and effTransProb is called)
  # it's much easier to do right here.
  # This kind of makes sense anyways since these transitions are exact only for ML inference
  # (ie they're effective transitions rather than true parts of the model).
  print "/*\n * transition_scores initialization: effective direct transitions \n * due to half-empty bifurcations \n */\n";
  foreach my $srcKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm}) {
    my $s = $self->cState ($srcKey);
    if (!$self->isNull ($s)) { next; }

    foreach my $destKey (sort {$self->stateCmp ($a, $b)} keys %{$self->tm->{$srcKey}}) {
      my $d = $self->cState ($destKey);
      if (!$self->isBifurc ($d)) { next; }
      
      my ($leftKey) = keys %{$self->bifurcMap->{$destKey}->{left}};
      my ($centerKey) = keys %{$self->bifurcMap->{$destKey}->{center}};

      # empty left child
      my $leftEmpty = 0;
      foreach my $leftChildKey (keys %{$self->tm->{$leftKey}}) {
	if ($self->isEnd ($self->cState ($leftChildKey))) { $leftEmpty = 1; last; }
      }
      if ($leftEmpty && ($self->isNull ($self->cState ($centerKey)))) {
	print "/*\n";
	print " * $srcKey -> ", Util->toAlphanumeric ($destKey), " -> (", Util->toAlphanumeric ($leftKey), ", ", Util->toAlphanumeric ($centerKey), ")\n";
	print " *   ", Util->toAlphanumeric ($leftKey), " -> End\n";
	print " *   ", Util->toAlphanumeric ($centerKey), "\n";
	print " * giving an effective transition\n";
	print " *   $srcKey -> ", Util->toAlphanumeric ($centerKey), "\n";
	print " */\n";

	print "transition_scores.transition ($srcKey, ", Util->toAlphanumeric ($centerKey), ") = ";
	print "ScorePMul (transition_scores.transition ($srcKey, ", Util->toAlphanumeric ($destKey), "), transition_scores.end\[", Util->toAlphanumeric ($leftKey), "\]);\n";
	print "\n";
      }

      # empty center child
      my $centerEmpty = 0;
      foreach my $centerChildKey (keys %{$self->tm->{$centerKey}}) {
	if ($self->isEnd ($self->cState ($centerChildKey))) { $centerEmpty = 1; last; }
      }
      if ($centerEmpty && ($self->isNull ($self->cState ($leftKey)))) {
	print "/* \n";
	print " * $srcKey -> ", Util->toAlphanumeric ($destKey), " -> (", Util->toAlphanumeric ($centerKey), ", ", Util->toAlphanumeric ($leftKey), ")\n";
	print " *   ", Util->toAlphanumeric ($centerKey), " -> End\n";
	print " *   ", Util->toAlphanumeric ($leftKey), "\n";
	print " * giving an effective transition\n";
	print " *   $srcKey -> ", Util->toAlphanumeric ($leftKey), "\n";
	print " */\n";

	print "transition_scores.transition ($srcKey, ", Util->toAlphanumeric ($leftKey), ") = ";
	print "ScorePMul (transition_scores.transition ($srcKey, ", Util->toAlphanumeric ($destKey), "), transition_scores.end\[", Util->toAlphanumeric ($centerKey), "\]);\n";
	print "\n";
      }

    }
  
  }
  
}

sub show {
  my ($self) = @_;

  # write formatted TripletSCFG if requested
  if ($self->directives_write) { $self->writeFormatted(); }

  # redundant states
  if ($self->directives_redundant) { $self->showRedundant ($self->tm); }

}

1
