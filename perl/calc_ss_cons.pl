#!/usr/bin/perl -w

use Stockholm::Database;

# command-line params
my $liberal = 0;
my $drop = 0;

# usage
my $usage = "Usage: $0 <Stockholm database>
            [-h, --help] Print this message.
            [-l, --liberal] Take the union of all #=GR SS lines to calculate consensus SS.
                         Conflicting base-pairs not annotated in SS_cons line.
                         Default is conservative, ie to take the intersection.
            [-d, --drop] Drop columns with conflicting base-pair annotations from the alignment.

Use #=GR SS lines to calculate a consensus secondary structure and add a #=GC SS_cons line.
By default a base-pair is annotated in the consensus structure only if all #=GR lines agree.
#=GR SS lines are unmodified.
\n";

my @argv;
while (@ARGV) {
  my $arg = shift;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif (($arg eq "-l") || ($arg eq "--liberal")) { $liberal = 1; }
    elsif (($arg eq "-d") || ($arg eq "--drop")) { $drop = 1; }
    else { die $usage; }
  } else {
    push @argv, $arg;
  }
}

unless (@argv) {
  @argv = ('-');
  warn "[waiting for alignments on standard input]\n";
}

my ($file) = @argv;
my $db = Stockholm::Database->from_file ($file);

my $cnt = 0;
for my $stock (@$db) {

  warn "[processing ",++$cnt,"/",scalar @$db," alignments]\n";
  unless ($stock->is_flush()) { die "File not flush.\n"; }

  my $cols = $stock->columns;

  my %bp = map (($_=>{}), @{$stock->seqname});
  for my $seqname (@{$stock->seqname}) {
    if (exists $stock->gr_SS->{$seqname}) {
      my $ss = $stock->gr_SS->{$seqname};
      $bp{$seqname} = parse_ss ($ss);
    } else { warn "[no SS line for $seqname]\n"; }
  }

  # take union of base-pairs
  my %bp_union;
  my %bp_union_r; # right-to-left lookup
  foreach my $seqbp (values %bp) {
    while (my ($lcol,$rcol) = each %$seqbp) {
      $bp_union{$lcol}->{$rcol} = 1; # allow for potential conflicts
      $bp_union_r{$rcol}->{$lcol} = 1;
    }
  }

  # liberal; look for non-conflicting base-pairs
  # (go through this even for conservative in case we want to drop conflicting base-pairs)
  my $ss_cons = "." x $cols;
  my $ss_tmp = "." x $cols;

  foreach my $lcol (keys %bp_union) {
    if (defined $bp_union_r{$lcol}) { next; } # if lcol annotated as > in some seq
    
    my @rcols = keys %{$bp_union{$lcol}};
    if (@rcols > 1) { next; } # if lcol has conflicting pairings with >
    
    else {
      my ($rcol) = @rcols;
      
      if (defined $bp_union{$rcol}) { next; } # if rcol annotated as < in some seq
      
      my @lcols = keys %{$bp_union_r{$rcol}};
      if (@lcols > 1) { next; } # if rcol has conflicting pairings with <
      elsif ($lcol == $lcols[0]) { # if unique left-right pairing
	if ((substr ($ss_tmp, $lcol, 1) ne ".") || (substr ($ss_tmp, $rcol, 1) ne ".")) { die "Oops, we've been here before: ($lcol $rcol)\n"; }
	substr ($ss_tmp, $lcol, 1) = "<";
	substr ($ss_tmp, $rcol, 1) = ">";
      } else { die "Unreachable\n"; }
    }
  }
  # now take care of case
  # X   < . > .
  # Y   . < . >
  # (so far it looks like '<<>>', but it really should be '....')
  # first parse the ss_tmp which we have now, then drop all basepairs which aren't in bp_union
  my $parse = parse_ss ($ss_tmp);
  while (my ($lcol,$rcol) = each %$parse) {
    if (!defined $bp_union{$lcol}->{$rcol}) { next; } # false base-pair
    else {
      substr ($ss_cons, $lcol, 1) = "<";
      substr ($ss_cons, $rcol, 1) = ">";
    }
  }

  # now calculate columns with conflicting base-pairs in case we want to drop them
  # conflicting base-pairs are those which don't appear in liberal SS_cons but do in a #=GR SS line
  my @conflicts;
  for (my $col = 0; $col < $cols; ++$col) {
    my $c = substr ($ss_cons, $col, 1);
    if (is_lchar ($c) || is_rchar ($c)) { next; }
    else {
      foreach my $seqname (keys %{$stock->gr_SS}) {
	my $ss = $stock->gr_SS->{$seqname};
	if (is_lchar (substr ($ss, $col, 1)) || is_rchar (substr ($ss, $col, 1))) {
	  push @conflicts, $col; last;
	}
      }
    }
  }

  # if conservative, then re-calculate $ss_cons as appropriate
  my %revbp = map (($bp{$_}=>$_), keys %bp);
  if (!$liberal) {
    $ss_cons = "." x $cols;
    foreach my $lcol (keys %bp_union) {
      foreach my $rcol (keys %{$bp_union{$lcol}}) {
	my $unanimous = 1;
	for my $seqbp (values %bp) {
	  if (!exists $seqbp->{$lcol} || $seqbp->{$lcol} != $rcol) {
	    $unanimous = 0;
	    last;
	  }
	}
	if ($unanimous || $liberal) {
	  substr ($ss_cons, $lcol, 1) = "<";
	  substr ($ss_cons, $rcol, 1) = ">";
	}
      }
    }
  }

  # store SS_cons
  assert_ss_valid ($stock->gc_SS_cons);
  $stock->gc_SS_cons ($ss_cons);

  # drop conflicting cols if requested
  if ($drop) { $stock->drop_columns (@conflicts); }
  assert_ss_valid ($stock->gc_SS_cons);
  while (my ($seqname, $ss) = each %{$stock->gr_SS}) {
    assert_ss_valid ($ss);
  }

  print $stock->to_string;
}

sub is_lchar {
  my ($c) = @_;
  return $c eq '<' || $c eq '[' || $c eq '(';
}

sub is_rchar {
  my ($c) = @_;
  return $c eq '>' || $c eq ']' || $c eq ')';
}

# Parse secondary structure.  Return hash mapping lcol -> rcol.
# Die if invalid ss.
sub parse_ss {
  my ($ss) = @_;

  my %bp;
  my @lcol;
  for (my $col = 0; $col < length $ss; ++$col) {
    my $c = substr ($ss, $col, 1);
    if (is_lchar ($c)) {
      push @lcol, $col;
    } elsif (is_rchar ($c)) {
      my $lcol = pop @lcol;
      if (!defined $lcol) { die "Invalid secondary structure:\n",$ss,"\n"; }
      $bp{$lcol} = $col;
    }
  }
  if (@lcol) { die "Invalid secondary structure:\n",$ss,"\n"; }

  return \%bp;
}

# Die if invalid ss.
sub assert_ss_valid {
  my ($ss) = @_;

  my @lcol;
  for (my $col = 0; $col < length $ss; ++$col) {
    my $c = substr ($ss, $col, 1);
    if (is_lchar ($c)) {
      push @lcol, $col;
    } elsif (is_rchar ($c)) {
      my $lcol = pop @lcol;
      if (!defined $lcol) { die "Invalid secondary structure:\n",$ss,"\n"; }
    }
  }
  if (@lcol) { die "Invalid secondary structure:\n",$ss,"\n"; }

  return 1;
}
