#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm::Database;

use strict;

# Stockholm params
my $SS = "SS";
my $SS_cons = "SS_cons";

# cmd-line params
my $dobp = 1;
my $nogaps = 0;
my $strictstems = 0;
my $allcons = 0;
my $looplen = 2;
my $mini = 0;
my $frac = 1;

# usage
my $usage = "Usage: $0 <Stockholm database>
            [-h, --help] Print this message.
            [--nogaps] Skip columns with any gaps.
            [--strictstems] Enforce Xrate-style 'strict-gaps' for base-pairs (no characters paired with gaps);
                     disallow base-pairs aligned with single-stranded characters.
            [--singlestranded] Extract single-stranded columns (default is to extract base-paired columns).
            [--allcons] Print all columns of alignment with conserved secondary structure (enforces strictstems).
            [--mini] Print miniature alignments (1 per stem or 1 per 10 single-stranded nt).
                     Ignores stems of length 1.
            [--looplen <int>] Create closing loops of the specified length (default is 2) by sampling from single-stranded columns.
                     Only effective in combination with --mini.
            [--frac <fraction>] Only process a random fraction of the alignments in the database.

Extracts base-paired or single-stranded columns from the input alignment based on the SS_cons line.
Skips all-gap columns.
Base-paired columns are NOT re-arranged (so e.g. stacking is preserved); single-stranded columns are
rearranged only if the --looplen option is used.
Single-stranded columns are extracted only if they aren't annotated as base-paired in any #=GR SS line as well.
\n";

my @argv;
while (@ARGV) {
  my $arg = shift;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif (($arg eq "--nogaps")) { $nogaps = 1; }
    elsif (($arg eq "--strictstems")) { $strictstems = 1; }
    elsif (($arg eq "--singlestranded")) { $dobp = 0; }
    elsif (($arg eq "--allcons")) { $allcons = 1; $strictstems = 1; }
    elsif (($arg eq "--mini")) { $mini = 1; }
    elsif (($arg eq "--looplen")) { $looplen = shift; }
    elsif (($arg eq "--frac")) { $frac = shift; }
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

  unless (rand(1) < $frac) {
    warn "[skipping ",++$cnt,"/",scalar @$db," alignments]\n";
    next;
  }
  warn "[processing ",++$cnt,"/",scalar @$db," alignments]\n";

  my $cols = $stock->columns;
  my $ss_cons = $stock->gc_($SS_cons);
  unless ($stock->is_flush()) { die "File not flush.  Is the SS_cons line missing?\n"; }
  assert_ss_valid ($ss_cons);
  while (my ($seqname, $ss) = each %{$stock->gr_($SS)}) {
    assert_ss_valid ($ss);
  }

  # @bp are base-paired in SS_cons (array refs)
  # @bp_discarded are base-paired in SS_cons but violate strict stem gaps
  # @ss are single-stranded in both #=GC SS_cons and #=GR SS lines
  # @ss_discarded are single-stranded in #=GC SS_cons but not in #=GR SS lines
  # @stem is array of arrayrefs giving indices of each stem in @bp
  my (@lpos, @bp, @bp_discarded, @ss, @ss_discarded, @stem);
  for (my $col = 0; $col < $cols; ++$col) {
    my $ss_char = substr ($ss_cons, $col, 1);
    if (is_lchar ($ss_char)) {
      push @lpos, $col;
    } elsif (is_rchar ($ss_char)) {
      my $l = pop @lpos;
      if (!defined $l) { die "Stack empty; invalid structure!\n"; }
      # enforce strict stem gaps if requested
      if ($strictstems && violates_strict_pairing ($stock, $l, $col)) { push @bp_discarded, [$l,$col]; next; }
      if ($nogaps && (has_gaps_col ($stock->get_column ($l)) || has_gaps_col ($stock->get_column ($col)))) { push @bp_discarded, [$l,$col]; next; } # don't store columns with gaps if so requested

      push @bp, [$l, $col];

      # keep track of stems
      if (is_lchar (substr ($ss_cons, $l+1, 1)) && is_rchar (substr ($ss_cons, $col-1, 1))) {
	# ...stacked...
	my $bp_index = @bp - 1; # index of current base-pair within @bp
	unless (@stem && $stem[-1]->[0] == $bp_index - 1) {  # was previous bp part of a stack? ($stem[-1] gets last element of @stem)
	  push @stem, [$bp_index - 1];  # if not, start a new stack
	}
	unshift @{$stem[-1]}, $bp_index; # store index of current base-pair 
      }      
    } else { # single-stranded

      # check that there's no per-seq base-pair annotation
      my $found_gr_ss = 0;
      foreach my $seqname (keys %{$stock->gr_($SS)}) {
	my $ss = $stock->gr_($SS)->{$seqname};
	if (is_lchar (substr ($ss, $col, 1)) || is_rchar (substr ($ss, $col, 1))) {
	  $found_gr_ss = 1;
	  last;
	}
      }
      if ($found_gr_ss) {
	push @ss_discarded, $col;
      } else {
	if (is_allgaps_col ($stock->get_column ($col))) { push @ss_discarded, $col; next; } # don't store all-gap columns
	if ($nogaps && has_gaps_col ($stock->get_column ($col))) { push @ss_discarded, $col; next; } # don't store columns with gaps if so requested
	push @ss, $col;
      }
    }
  }

  # sanity check
  die "Oops" unless (2*scalar @bp) + (2*scalar @bp_discarded) + @ss + @ss_discarded == $cols;

  # print miniature alignments if requested
  if ($mini) {

    # base-paired
    if ($dobp) {
      foreach my $stem (@stem) { 

	# sample from @ss to create loop sequence for each stem
	my @bp_slice = @bp[@$stem];
	my @loop_slice;
	for (my $l = 0; $l < $looplen; ++$l) {
	  push @loop_slice, @ss[int(rand(@ss))];
	}

	my @stem_slice;
	push @stem_slice, map ($_->[0], @bp_slice);
	push @stem_slice, @loop_slice;
	push @stem_slice, map ($_->[1], reverse (@bp_slice));
		    
	print_sliced_copy ($stock, \@stem_slice);
      }
    }

    # single-stranded
    else {
      my $i = 0;
      for ($i = 0; $i < @ss-10; $i += 10) { # print by 10 nt
	my @mini_slice = @ss[$i..$i+10];
	print_sliced_copy ($stock, \@mini_slice);
      }
      my @mini_slice = @ss[$i+1..scalar(@ss)-1]; # catch the last bits
      print_sliced_copy ($stock, \@mini_slice);
    }

  }

  # else just print the regular stuff
  else {

    my @slice_cols;
    if ($allcons) {
      push @slice_cols, sort { $a <=> $b } (map (@$_, @bp), @ss);
    } else {
      if ($dobp) {
	push @slice_cols, sort { $a <=> $b } (map (@$_, @bp));
      } else {
	push @slice_cols, sort { $a <=> $b } (@ss);
      }
    }
    print_sliced ($stock, \@slice_cols);

  }

}

# is a character a gap?
sub is_gap {
  my ($c) = @_;
  my $gapChars = '-._';
  return $c =~ /[$gapChars]/;
}

# does a column have any gaps?
sub has_gaps_col {
  my ($col) = @_;
  for (my $r = 0; $r < length $col; ++$r) { 
    if (is_gap (substr ($col,$r,1))) { return 1; }
  }
  return 0;
}

# is a column all gaps?
sub is_allgaps_col {
  my ($col) = @_;
  for (my $r = 0; $r < length $col; ++$r) { 
    if (!is_gap (substr ($col,$r,1))) { return 0; }
  }
  return 1;
}

# Do 2 columns annotated as base-paired in $SS_cons:
# 1) violate Xrate-style strict-gaps (no characters paired w/ gaps)?
# 2) have characters annotated as single-stranded in #=GR $SS lines?
# NB 0-based column indexing.
sub violates_strict_pairing {
  my ($stock, $lcol, $rcol) = @_;
  
  foreach my $seqname (keys %{$stock->gr_($SS)}) {
    my $seq = $stock->seqdata->{$seqname};
    my ($lc, $rc) = (substr ($seq, $lcol, 1), substr ($seq, $rcol, 1));
    my $ss = $stock->gr_($SS)->{$seqname};
    # check strict-gaps: no characters paired with gaps
    if ((is_gap ($lc) && !is_gap ($rc)) || (!is_gap ($lc) && is_gap ($rc))) { return 1; }
    # check for no single-stranded annotations for non-gap characters
    if (!is_gap ($lc) && !is_lchar (substr ($ss, $lcol, 1))) { return 1; }
    if (!is_gap ($rc) && !is_rchar (substr ($ss, $rcol, 1))) { return 1; }
  }
  
  return 0;
}

sub is_lchar {
  my ($c) = @_;
  return $c eq '<' || $c eq '[';
}

sub is_rchar {
  my ($c) = @_;
  return $c eq '>' || $c eq ']';
}

sub slice_seq {
  my ($textRef, $slice) = @_;
  $$textRef = join ("", (split (//, $$textRef)) [@$slice]);
}

sub print_sliced {
  my ($stock, $slice) = @_;

  # slice alignment columns
  for my $seqname (@{$stock->seqname}) {
    slice_seq (\$stock->seqdata->{$seqname}, $slice);
    while (my ($tag, $hash) = each %{$stock->gr}) {
      if (exists ($hash->{$seqname})) {
	slice_seq (\$hash->{$seqname}, $slice);
      }
    }
  }
  for my $tag (keys %{$stock->gc}) {
    slice_seq (\$stock->gc->{$tag}, $slice);
  }

  # sanity check that consensus structure and individual structures are valid
  assert_ss_valid ($stock->gc_($SS_cons));
  while (my ($seqname, $ss) = each %{$stock->gr_($SS)}) {
    assert_ss_valid ($ss);
  }

  # print
  print $stock->to_string;
  
}

# like print_sliced, but copy to avoid modifying the original alignment
sub print_sliced_copy {
  my ($stock, $slice) = @_;
  my $copy = $stock->copy;
  print_sliced ($copy, $slice);
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
