#!/usr/bin/perl -w

# NB: The End state /must/ be named 'e' in all .tt files.

package StateTyping;

use strict;
use vars '@ISA';

# Is a state type of type start?
sub isStart {
    my ($class, $t) = @_;
    return ($t eq 's');
}

# Is a state type of type end?
sub isEnd {
    my ($class, $t) = @_;
    return ($t eq 'e');
}

# Is a state type of type wait?
sub isWait {
    my ($class, $t) = @_;
    return ($t eq 'w');
}

# Is a state type of type insert?
sub isInsert {
    my ($class, $t) = @_;
    return ($t eq 'i');
}

# Is a state type of type match?
sub isMatch {
    my ($class, $t) = @_;
    return ($t eq 'm');
}

# Is a state type of type bifurcation?
sub isBifurc {
    my ($class, $t) = @_;
    return ($class->isBifurcInsert ($t) || $class->isBifurcMatch ($t));
}

sub isBifurcInsert {
  my ($class, $t) = @_;
  return ($t eq 'bi');
}

sub isBifurcMatch {
  my ($class, $t) = @_;
  return ($t eq 'bm');
}

# Return the (unique) End state.
sub endState {
  my ($class) = @_;
  return 'e';
}


########################################################################
package TreeTransducer;

use Carp;

use strict;
use vars '@ISA';

sub new {
  my ($class) = @_;
  my $self = {
	      'stateTyping'  => {},
	      'absorbProfiling' => {},
	      'emitProfiling' => {},
	      'emitDist' => {},
	      'tm'      => {},
	      'effTM'      => {},
	      'startState' => "", # the "overall" start state
	      'bifurcMap' => {},
	      'stateSorting'       => {},
	     };
  bless $self, $class;

  return $self;
}

# catch methods
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD; # $AUTOLOAD contains the fully qualified name of the original subroutine
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq 'DESTROY';

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
  
  $self->effTM ($self->buildEffTM());
}

# Various checks that the transducer is valid.
# Some checks are relevant only for branch machines.
sub _checkValid {
  my ($self) = @_;

  # Check that all states of type Insert, BifurcInsert and Match, BifurcMatch have
  # absorb and emit profiles as appropriate.
  foreach my $s (@{$self->states()}) {
    if ($self->isInsert ($s)) {
      my $p = $self->emitProfile ($s);
      if (!$p || !($p eq 'l' || $p eq 'r' || $p eq 'lr')) { # insert states must have emit profiles
	croak ("Empty or invalid emit profile found for state '$s'.\n");
      }
      if ($p && !exists $self->emitDist->{$s}) {
	croak ("No emit distribution defined for state '$s'.\n");
      }
    }
    elsif ($self->isMatch ($s)) {
      my $a = $self->absorbProfile ($s);
      my $p = $self->emitProfile ($s);
      if (!$a || !($a eq 'l' || $a eq 'r' || $a eq 'lr')) { # match states must have absorb profiles
	croak ("Empty or invalid absorb profile found for state '$s'.\n");
      }
      if ($p && !($p eq 'l' || $p eq 'r' || $p eq 'lr')) { # match states can have empty emit profiles (e.g. delete states)
	croak ("Invalid emit profile found for state '$s'.\n");
      }
      if ($p && !exists $self->emitDist->{$s}) {
	croak ("No emit distribution defined for state '$s'.\n");
      }
    }
    elsif ($self->isBifurc ($s)) {
      if (!exists $self->emitProfiling->{$s}) {
	croak ("Empty emit profile found for Bifurc state '$s'.\n");
      }
      if ($self->isBifurcInsert ($s)) {
	if (exists $self->absorbProfiling->{$s}) {
	  croak ("BifurcInsert state '$s' cannot have an absorb profile.\n");
	}
      } elsif ($self->isBifurcMatch ($s)) {
	if (!exists $self->absorbProfiling->{$s}) {
	  croak ("BifurcMatch state '$s' must have a non-empty absorb profile.\n");
	}
      }
    }
  }

  # Checks on the transition matrix and bifurcations:
  #  - Every dest state in the TM must either
  #     1) be a source state (if not a bifurc)
  #     2) have children (if a bifurc)
  #  - Bifurcation states can't appear in the TM as source states.
  #  - The children of bifurc states must:
  #     1) be source states in the TM.
  #     2) be of type Start or End
  foreach my $src (@{$self->states()}) {
    if ($self->isBifurc ($src)) {
      if (exists $self->tm->{$src}) {
	croak ("bifurcation state '$src' shouldn't be present in the TM as a source state.\n");
      }
      if (!ref $self->bifurcChildren ($src)) {
	croak ("bifurcation state '$src' has no children.\n");
      }
      foreach my $dest (@{$self->destStates ($src)}) {
	if (!($self->isStart ($dest) || $self->isEnd ($dest))) {
	  croak ("bifurcation state '$src' has a child '$dest' which isn't of type Start or End.\n");
	}
      }
    }
    foreach my $dest (@{$self->destStates ($src)}) {
      if ($self->isBifurc ($dest) && !ref $self->bifurcChildren ($dest)) {
	croak ("Invalid: bifurcation state '$dest' has no children.\n");
      }
      if (!$self->isBifurc ($dest) && !ref($self->tm->{$dest})) {
	croak ("Invalid TM: exists transition $src -> $dest, but '$dest' isn't present as a source state.\n");
      }
    }
  }


  # check that each state has a defined sorting parameter
  foreach my $s (@{$self->states()}) {
    if (!exists $self->stateSorting->{$s} || !($self->stateSorting->{$s} =~ /\d+/)) {
      croak ("No sorting parameter found for state '$s'.\n");
    }
  }


}

sub stateType {
  my ($self, $state) = @_;
  return $self->stateTyping->{$state};
}

sub isStart {
  my ($self, $state) = @_;
  return (StateTyping->isStart ($self->stateType ($state)));
}

sub isEnd {
  my ($self, $state) = @_;
  return (StateTyping->isEnd ($self->stateType ($state)));
}

sub isWait {
  my ($self, $state) = @_;
  return (StateTyping->isWait ($self->stateType ($state)));
}

sub isInsert {
  my ($self, $state) = @_;
  return (StateTyping->isInsert ($self->stateType ($state)));
}

sub isMatch {
  my ($self, $state) = @_;
  return (StateTyping->isMatch ($self->stateType ($state)));
}

sub isBifurc {
  my ($self, $state) = @_;
  return (StateTyping->isBifurc ($self->stateType ($state)));
}

sub isBifurcInsert {
  my ($self, $state) = @_;
  return (StateTyping->isBifurcInsert ($self->stateType ($state)));
}

sub isBifurcMatch {
  my ($self, $state) = @_;
  return (StateTyping->isBifurcMatch ($self->stateType ($state)));
}

sub endState {
  my ($class) = @_;
  return StateTyping->endState();
}

# Given a state of type Match or Bifurc, get
# the corresponding absorption profile.
sub absorbProfile {
  my ($self, $s) = @_;
  if ($self->isMatch ($s) || $self->isBifurcMatch ($s)) {
    return $self->absorbProfiling->{$s};
  }    
  warn "state '$s' is not of type Match or BifurcMatch.\n";
}

# Given a state of type Insert or Bifurc, get the corresponding emission profile.
sub emitProfile {
  my ($self, $s) = @_;
  if ($self->isInsert ($s) || $self->isMatch ($s) || $self->isBifurc ($s)) {
    return $self->emitProfiling->{$s};
  }
  warn "State '$s' is not of type Insert or BifurcInsert.\n";
}

# Given a source and destination state, get the corresponding transition probability
# from the tm (returns a string).  If a time parameter is passed, incorporate it using formatProb().
# Handles bifurcations too.
# This really should only be called internally to the class; generally we'll use
# effTransProb().
sub transProb {
  my ($self, $src, $dest, $time) = @_;

  # handle bifurcation case
  if ($self->isBifurc ($src)) {
    # check that it's a valid bifurcation
    my $ok = 0;
    foreach my $child (@{$self->destStates ($src)}) {
      if ($dest eq $child) { $ok = 1; }
    }
    unless ($ok) { croak ("Invalid bifurcation: '$src' has no child '$dest'.\n"); }
    return '(1)';
  }

  my $prob;
  if (exists $self->tm->{$src} and
      exists $self->tm->{$src}->{$dest}) {
    $prob = $self->tm->{$src}->{$dest};
    if (defined $time) { $prob = $self->formatProb ($prob, $time); }
  } else {
    return '(0)';
  }
  $prob = '(' . $prob . ')';

  return $prob;
}


# Returns the effective transition probability between 2 states.
# See buildEffTM for details.
sub effTransProb {
  my ($self, $src, $dest, $time) = @_;

  # handle bifurcation case
  if ($self->isBifurc ($src)) {
    # check that it's a valid bifurcation
    my $ok = 0;
    foreach my $child (@{$self->destStates ($src)}) {
      if ($dest eq $child) { $ok = 1; }
    }
    unless ($ok) { croak ("Invalid bifurcation: $src has no child $dest.\n"); }
    return '(1)';
  }

  my $prob;
  if (exists $self->effTM->{$src} and
      exists $self->effTM->{$src}->{$dest}) {
    $prob = $self->effTM->{$src}->{$dest};
    if (defined $time) { $prob = $self->formatProb ($prob, $time); }
  } else {
    return '(0)';
  }

  return $prob;
}

# For each pair of states (dest state not a Wait state),
# store the effective transition probability.
# This is defined as:
# 1) the prob as defined in the original tm
# 2) the path probability between 'src' and 'dest' where the path 
#      consists of bifurcation states which have only 1 non-End child
#      (e.g.  l -> B -> (s e e) replaced by an effective l -> s)
#      NB: 'dest' must be of type Start
#      These bifurcations are later "pruned" out of the composed machine.
# 3) the sum over all paths probabilities between 'src' and 'dest'
#      where the paths consist only of single Wait states
#      (relevant only for Branch machines)
# Wraps everything nicely in parentheses.
sub buildEffTM {
  my ($self) = @_;

  my $effTM = {};

  # First construct the effTM considering only 1) and 2).
  foreach my $src (@{$self->states()}) {
    # Bifurc states shouldn't be present as source states in the tm
    if ($self->isBifurc ($src)) { next; }

    foreach my $dest (@{$self->destStates ($src)}) {
      # 1) if a direct transition exists, store it
      $effTM->{$src}->{$dest} = '(' . $self->tm->{$src}->{$dest} . ')';

      # 2) examine dest bifurc states to look for ones with 2 End children      
      unless ($self->isBifurc ($dest)) { next; }
      my $cntEnd = 0;
      my $child;
      foreach my $c (@{$self->bifurcChildren ($dest)}) {
	if ($self->isEnd ($c)) { $cntEnd += 1; }
	else { $child = $c; }
      }
      if ($cntEnd == 2) {
	$effTM->{$src}->{$child} = '(' . $self->transProb ($src, $dest) . '*' . $self->transProb ($dest, $child) . ')';
      }
    }
  }

  # Now do the summation over paths for 3), being careful to use the transitions and transition probabilities just defined in $effTM
  foreach my $src (@{$self->states()}) {
    # Bifurc states shouldn't be present as source states in the tm
    if ($self->isBifurc ($src)) { next; }

    foreach my $dest (@{$self->states()}) {
      # if this transition exists in the tm, continue
      if (exists $self->tm->{$src}->{$dest}) { next; }

      my @factors;
      foreach my $w (@{$self->destStates ($src)}) {
	unless ($self->isWait ($w)) { next; }
	if (exists $effTM->{$w} and # this first check shouldn't be necessary, but just being careful...
	    exists $effTM->{$w}->{$dest}) { # note here that we MUST refer to effTM, as it has the transitions due to pruned bifurcs as well
	  my $prob = '(' . $effTM->{$src}->{$w} . "*" . $effTM->{$w}->{$dest} . ')'; # here we MUST refer to the effTM probs rather than the tm ones!
	  push (@factors, $prob);
	}
      }
      if (scalar @factors) {
	if (scalar (@factors) == 1) { $effTM->{$src}->{$dest} = $factors[0]; }
	else { $effTM->{$src}->{$dest} = '(' . join (" + ", @factors) . ')'; }
      }
    }
  }

  return $effTM;
}

# Formats a transition probability string: replaces each alphanumeric 
# parameter 'a(t)' (the '(t)' means that it depends on time) in the transition
# probability with 'a(time)'.
sub formatProb {
  my ($self, $p, $time) = @_;

  if ($p =~ /(\w*[a-zA-Z]+\w*)\(t\)/) {
    $p =~ s/\(t\)/\($time\)/g;
  }

  return $p;
}

# Formats an emit dist string: e.g. replaces m_t with m_2.
sub formatEmitDist {
  my ($self, $d, $n) = @_;

  if ($d =~ /(\w*[a-zA-Z]+\w*)_t/) {
    $d =~ s/_t/_$n/g;
  }

  return $d;
}

# Given a bifurcation state, return a ref to an array of its child states.
sub bifurcChildren {
  my ($self, $src) = @_;
  if (!$self->isBifurc ($src)) { croak ("State '$src' isn't a bifurcation state.\n"); }
  return $self->bifurcMap->{$src};
}

# Get ref to array of all states.
sub states {
  my ($self) = @_;
  my @states = keys %{$self->tm};
  push (@states, @{$self->bifurcStates()});
  return \@states;
}

sub bifurcStates {
  my ($self) = @_;
  my @states = keys %{$self->bifurcMap};
  return \@states;
}

# Given a source state, return a ref to an array of
# the possible destination states.  If the source state
# is a bifurcation state, return a ref to an array of 
# the corresponding child states.
sub destStates {
  my ($self, $src) = @_;
  my @dest;

  # handle bifurcations separately
  if ($self->isBifurc ($src)) {
    return $self->bifurcChildren ($src);
  } else {
    @dest = keys %{$self->tm->{$src}};
    return \@dest;
  }
}


# Initializes a TreeTransducer from a flat file.
sub from_file {
  my ($class, $filename) = @_;

  my $self = $class->new();

  if (!(-e $filename)) { croak "'$filename' doesn't exist.\n"; }

  warn "Opening '$filename'.\n";
  open (FILE, "<$filename");
  warn "Parsing '$filename'.\n";
  my $label = undef;
  while (<FILE>) {
    # ignore comments and blank lines
    next if /^;/;
    next unless /\S+/;

    chomp();
    # if found label, record it
    if (/^\s*>\s*(\S+)/) {
      $label = $1;
      # print $label, "\n";
      next;
    }

    # process the line
    else {
      if (!defined $label) { warn "Ignoring: '$_'"; next; }

      # strip whitespace
      s/\s//g;

      if ($label =~ /statetyping/i) {
	if (/^(\S+)=(\w+)$/) {
	  my ($s, $t) = ($1, $2);
	  $self->stateTyping->{$s} = $t;
	} else { warn "Ignoring '$_'"; next; } }

      elsif ($label =~ /absorbprofiling/i) {
	if (/^(\S+)=(\S+)$/) {
	  my ($s, $p) = ($1, $2);
	  $self->absorbProfiling->{$s} = $p;
	} else { warn "Ignoring '$_'"; next; } }

      elsif ($label =~ /emitdist/i) {
	if (/^(\S+)=(\S+)$/) {
	  my ($s, $dist) = ($1, $2);
	  $self->emitDist->{$s} = $dist;
      }  else { warn "Ignoring '$_'"; next; } }

      elsif ($label =~ /emitprofiling/i) {
	if (/^(\S+)=(\S*)$/) {
	  my ($s, $p) = ($1, $2);
	  $self->emitProfiling->{$s} = $p;
	} else { warn "Ignoring '$_'"; next; } }

      elsif ($label =~ /tm/i) {
	if (/^(\S+)->(\S*)/) {
	  my ($from, $tolist) = ($1, $2);
	  my @grps = split (",", $tolist);
	  # if transitions defined
	  if (@grps) {
	    foreach my $grp (@grps) {
	      if ($grp =~ /^(\S+)=(\S+)$/) {
		my ($to, $prob) = ($1, $2);
		$self->tm->{$from}->{$to} = $prob;
	      } else { warn "Ignoring '$grp'"; next; }
	    }
	  }
	  # if not (end state)
	  else {
	    $self->tm->{$from} = {};
	  }
	  # store the "overall" start state
	  if (!$self->startState) {
	    if (!$self->isStart ($from)) { croak "'$from' is not a valid Start state.\n"; }
	    $self->startState ($from);
	  }
	} else { warn "Ignoring '$_'"; next; } }

      elsif ($label =~ /bifurc/i) {
	if (/^(\S+)->\((\S+)\)/) {
	  my ($b, $childlist) = ($1, $2);
	  my @children = split (",", $childlist);
	  if (scalar (@children) eq 3) {
	    push (@{$self->bifurcMap->{$b}}, @children);
	  } else { warn "3 child states must be specified (although empty states are ok for e.g. left bifurcations; ignoring '$_'"; next; }
	} else { warn "Ignorning '$_'"; next; } }

      elsif ($label =~ /statesorting/i) {
	if (/^(\S+)=(\d+\.?\d*)$/) {
	  my ($s, $scale) = ($1, $2);
	  $self->stateSorting->{$s} = $scale;
	} else { warn "Ignoring '$_'"; next; } }

      else { warn "I don't recognize the label '$label'.  Ignoring '$_'"; next; }
    }
  }
  close (FILE);

  $self->_initialize();
  $self->_checkValid();

  return $self;
}

sub showStatistics {
  my ($self) = @_;

  my $numStates = scalar (@{$self->states()});
  print ";   $numStates total states\n";
  print ";   ", scalar (@{$self->bifurcStates()}), " deterministic bifurcation states\n";
}

1
