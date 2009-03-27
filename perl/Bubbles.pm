#!/usr/bin/perl -w

package Bubbles;

use strict;
use vars '@ISA';

use Carp;

my @nuc =     qw/ X A C G T /;
my @nuc_rna = qw/ X A C G U /;
my @pur =     qw/ X 1 0 1 0 /; # purine is 1

my @dinuc =     qw/ X AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT /;
my @dinuc_rna = qw/ X AA AC AG AU CA CC CG CU GA GC GG GU UA UC UG UU /;
my @wc =        qw/ X 0  0  0  1  0  0  1  0  0  1  0  0  1  0  0  0  /; # Watson-Crick pairings
my @canonical = qw/ X 0  0  0  1  0  0  1  0  0  1  0  1  1  0  1  0  /; # canonical pairings

my @codon = qw/ X AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAC TAT TCA TCC TCG TCT TGC TGG TGT TTA TTC TTG TTT /;
my @codonAA = qw/ X K N K N T T T T R S R S I I M I Q H Q H P P P P R R R R L L L L E D E D A A A A G G G G V V V V Y Y S S S S C W C L F L F /;

my @aa = qw/ X A C D E F G H I K L M N P Q R S T V W Y /;
my @aaAromatic = qw/F W Y/; # hydrophobic as well
my @aaHydrophobic = qw/A C G I L M P V/; # non-aromatic
my @aaHydrophilic = qw/D E H K N Q R S T/; # includes charged aa 

sub new {
  my ($class, $directives, $N) = @_;
  my $self = {
	      'directives' => {}, # command-line flags (graphics directives); curly brackets mean hash reference
	      'N' => "",
	      'nt1' => [], # square brackets mean array reference
	      'nt2' => [],
	      'nt3' => [],
	      'ntcode1' => [],
	      'ntcode2' => [],
	      'ntcode3' => [],
	      'sorting' => [],
	      'grid' => [],
	     };
  bless $self, $class;
  $self->directives($directives); # set hashref
  $self->N($N);

  return $self->_initialize;
}

sub _initialize {
  my ($self) = @_;

  # stupid hack to override canonical base pair defs, if the user chooses to
  @wc = @canonical if ($self->directives->{'gupairs'});

  # choose ordering
  if    ($self->N == 4) {
    @nuc = @nuc_rna if ($self->directives->{'rnamode'});  # stupid hack to get correct alphabet
    $self->sorting(['X','1','2','3','4']); } #'X' is a dummy 0th element making indexing easier later on!
  elsif ($self->N == 16) {
    @dinuc = @dinuc_rna if ($self->directives->{'rnamode'});  # stupid hack to get correct alphabet

    # old way, where hardcoding was necessary [AVU]
    #$self->nt1(['X','A','A','A','A','C','C','C','C','G','G','G','G','T','T','T','T']);
    #$self->nt2(['X','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T']);

    # new way (less hardcoding, now you only have to change alphabet in the global declarations) [AVU]
    my @nt1 = ('X');  my @nt2 = ('X');   # use 'X' placeholder to allow 1-based indexing
    foreach my $pair (@dinuc[1..$#dinuc]) {
      push (@nt1, substr ($pair, 0, 1));
      push (@nt2, substr ($pair, 1, 1));
    }
    $self->nt1(\@nt1);
    $self->nt2(\@nt2);
    $self->ntcode1(['X','0','0','0','0','1','1','1','1','2','2','2','2','3','3','3','3']);
    $self->ntcode2(['X','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3']);
    $self->sorting(['X','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16']); }
  elsif (($self->N == 61) || ($self->N == 64)) {
    $self->nt1(['X','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','T','T','T','T','T','T','T','T','T','T','T','T','T']);
    $self->nt2(['X','A','A','A','A','C','C','C','C','G','G','G','G','T','T','T','T','A','A','A','A','C','C','C','C','G','G','G','G','T','T','T','T','A','A','A','A','C','C','C','C','G','G','G','G','T','T','T','T','A','A','C','C','C','C','G','G','G','T','T','T','T']);
    $self->nt3(['X','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','A','C','G','T','C','T','A','C','G','T','C','G','T','A','C','G','T']);
    $self->ntcode1(['X','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','3','3','3','3','3','3','3','3','3','3','3','3','3']);
    $self->ntcode2(['X','0','0','0','0','1','1','1','1','2','2','2','2','3','3','3','3','0','0','0','0','1','1','1','1','2','2','2','2','3','3','3','3','0','0','0','0','1','1','1','1','2','2','2','2','3','3','3','3','0','0','1','1','1','1','2','2','2','3','3','3','3']);
    $self->ntcode3(['X','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','0','1','2','3','1','3','0','1','2','3','1','2','3','0','1','2','3']);
    if (!$self->directives_resort) { # alphabetical sorting
      $self->sorting(['X','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61']); }
    elsif ($self->directives_resort == 1) { # genetic ordering
      $self->sorting(['X','40','39','41','38','36','35','37','34','44','43','45','42','32','31','33','30','24','23','25','22','20','19','21','18','28','27','29','26','16','15','17','14','56','55','57','54','52','51','53','50','60','59','61','58','48','47','49','46','10','9','7','6','8','5','12','13','11','3','2','4','1']); }
    elsif ($self->directives_resort == 2) { # used for re-ordering codons into synonymous blocks (based on alphabetical order of encoded amino acids)
      $self->sorting(['X','36','11','37','12','51','52','53','54','5','45','6','46','27','28','38','29','17','25','18','26','41','42','43','44','7','8','9','10','30','31','32','33','19','13','20','14','1','2','3','4','21','22','23','24','58','59','60','61','56','57','47','48','49','50','15','55','16','34','39','35','40']); }
    elsif ($self->directives_resort == 3) { # used for re-ordering codons into synonymous blocks (based on biochemical similarity)
      $self->sorting(['X','21','17','22','18','52','53','54','55','23','56','24','57','43','44','36','45','19','29','20','30','9','10','11','12','25','26','27','28','37','38','39','40','13','15','14','16','1','2','3','4','5','6','7','8','46','47','48','49','32','33','58','59','60','61','50','31','51','41','34','42','35']); }
    elsif ($self->directives_resort == 4) { # ordering due to Urbina and Tang (JME 2006) and Higgs
      $self->sorting(['X','41','40','42','39','27','26','28','25','56','55','57','54','11','10','12','9','37','36','38','35','23','22','24','21','52','51','53','50','7','6','8','5','45','44','46','43','31','30','32','29','60','59','61','58','15','14','16','13','34','33','19','18','20','17','48','49','47','3','2','4','1']); }
    else { croak "resort = ",$self->directives_resort," is not a valid flag.\n"; }
    if ($self->directives_newgrid) { # if grid by synonymous blocks
      if ($self->directives_resort == 2) { $self->grid(['4.5','10.5','12.5','14.5','16.5','18.5','20.5','24.5','26.5','29.5','35.5','37.5','38.5','40.5','44.5','50.5','54.5','55.5','57.5']); }
      elsif ($self->directives_resort == 3) { $self->grid(['4.5','8.5','12.5','14.5','16.5','18.5','20.5','22.5','28.5','30.5','31.5','33.5','35.5','36.5','42.5','45.5','49.5','51.5','55.5']); }
      else { $self->grid(['4.5','8.5','12.5','16.5','20.5','24.5','28.5','32.5','36.5','40.5','44.5','48.5','50.5','54.5','57.5']); } } }
  elsif ($self->N == 20) {
    $self->sorting(['X','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']); }
  else { croak "Huh?  Not a chain over nuc, dinuc, codons or a.a.!\n"; }

  return $self;
}


# catch methods
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD; # $AUTOLOAD contains the fully qualified name of the original subroutine
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq 'DESTROY';

  # check for directives accessor, e.g. $self->directives_rscale or $self->directives_('rscale')
  if ($sub =~ /^directives_(\S*)$/i) {
    my $flag = lc($1);
    $flag = shift @args unless (length $flag); # this catches the argument 'Ka' in the second example usage given above
    if (!defined $self->{"directives"}->{$flag}) {
      $self->{"directives"}->{$flag} = "";     # if no such flag exists, create one with the value ""
    }                                          # we therefore have to test 'if ($self->directives_newgrid)' rather than 'if (defined $self->directives_newgrid)'
    return $self->{"directives"}->{$flag};     # the second will always be true because AUTOLOAD will create an empty flag for us
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


## Draw inital probability distribution over states ##
sub drawInitDist {
  my ($self, $mp, $dat) = @_;
  my $N = $self->N;
  my ($xoffset, $yoffset) = ($self->directives_xoffset,$self->directives_yoffset);
  my $penscale = $self->directives_penscale;

  #draw grid (currently $newgrid stuff not implemented here)
  print $mp "% draw grid\n";
  print $mp "draw ($xoffset*u,($yoffset+1)*u)--(($N+1+$xoffset)*u,($yoffset+1)*u) withpen pencircle scaled $penscale withcolor 0.7white;\n"; # horizontal line
  print $mp "for i=1 upto $N:\n  draw (($xoffset+i)*u,$yoffset*u)--(($xoffset+i)*u,($yoffset+2)*u) withpen pencircle scaled $penscale withcolor 0.6white;\nendfor\n"; # vertical lines

  my $count = 0;
  my (@x, @y, @z);
  #for each row of input file...
  while (<$dat>) {
    #pick out data
    my @no = split;
    my ($state, $value) = ($no[0], $no[1]);
    #if a valid matrix element...
    unless (($state < 1) || ($state > $N)) {
      my ($col, $row) = ($self->sorting->[$state], $N);
      $col += $xoffset; # move to right if necessary to be flush with rate plot
      $row -= $yoffset; # minus because metapost definitions are such that y-axis is flipped (origin is top-left corner)

      #output most of the circle command, except the color information
      my $scaledvalue = 3 * $self->directives_pscale * sqrt($value);
      if ($self->directives_smooth) { # MetaPost seems to barf when $scaledvalue is in scientific notation; ugly printf hack to fix this
	print $mp "dullfullcircle(($col*u,$row*u),";
	printf $mp "%f*u,", $scaledvalue;
      }
      else {
	croak "prettyfullcircle stuff not implemented yet!\n";
#	print $mp "prettyfullcircle(($col*u,$row*u),";
#	printf $mp "%f*u,", $scaledvalue;
#	print $mp "$theta,$fraction,$smooth,";
      }
      #if needed, store values for re-drawing circles on top
      if ($self->directives_forceoutline) { $x[$count] = $col; $y[$count] = $row; $z[$count] = $scaledvalue; $count++; }

      if    ($N == 4) { print $mp "(0,0,0),$N);\n"; }
      elsif ($N == 16) { if ($wc[$state]) { print $mp "(0,1,0),$N);\n"; } # Watson-Crick paired = green
			 else { print $mp "(0,0,0),$N);\n"; } } # non-WC = black
      elsif (($N == 61) || ($N == 64)) { print $mp "(0,0,1),$N);\n";  } # all codons same color for now
      elsif ($N == 20) { 
	my $aminoAcid = $aa[$state];
	### scheme 0
	if ($self->directives_colsch == 0) {
	  if ( grep($aminoAcid eq $_, @aaAromatic) ) {
	    print $mp "(1,0,0),$N);\n"; # red
	  }
	  elsif ( grep($aminoAcid eq $_, @aaHydrophobic) ) {
	    print $mp "(0,1,0),$N);\n"; # green
	  }
	  elsif ( grep($aminoAcid eq $_, @aaHydrophilic) ) {
	    print $mp "(0,0,1),$N);\n"; # blue
	  }
	}
      } 
    }
  }

  #if needed, re-draw circles on top
  if (($self->directives_colsch < 100) && ($self->directives_forceoutline)) {
    for (my $i = 0; $i < $count; $i++) {
      print $mp "dullemptycircle(($x[$i]*u,$y[$i]*u),$z[$i]*u,$N);\n";
    }
  }

  # label axes -- "states" label currently disabled for aesthetics
  print $mp "label.lft(\"prob.\",($xoffset*u,($yoffset+1)*u));\n";
  #  if    ($N == 4)  { print $mp "label.top(\"states\",(($N/2+$xoffset)*u,($yoffset+2.5)*u));"; }
  #  elsif ($N == 16) { print $mp "label.top(\"states\",(($N/2+$xoffset)*u,($yoffset+3.5)*u));"; }
  #  elsif (($N == 61) || ($N == 64)) { print $mp "label.top(\"states\",(($N/2+$xoffset)*u,($yoffset+4.5)*u));"; }
  #  elsif ($N == 20) { print $mp "label.top(\"states\",(($N/2+$xoffset)*u,($yoffset+2.5)*u));"; }

  # label states
  for (my $state = 1 ; $state < $N+1; $state++) {
    if ($state > 61) { warn "Skipping codon '$state' (stop codons disallowed)."; next; } # because we don't look at stop codons
    my $tmp = $self->sorting->[$state]; $tmp += $xoffset;
    if    ($N == 4)  { print $mp "label.top(\"",$nuc[$state],"\",(${tmp}*u,($yoffset+1.9)*u));\n"; }
    elsif ($N == 16) { print $mp "label.top(\"",$self->nt1->[$state],"\",(${tmp}*u,($yoffset+2.7)*u));\n";
                       print $mp "label.top(\"",$self->nt2->[$state],"\",(${tmp}*u,($yoffset+1.9)*u));\n"; }
    elsif (($N == 61) || ($N == 64)) {
		       print $mp "label.top(\"",$self->nt1->[$state],"\",(${tmp}*u,($yoffset+3.5)*u));\n";
		       print $mp "label.top(\"",$self->nt2->[$state],"\",(${tmp}*u,($yoffset+2.7)*u));\n";
		       print $mp "label.top(\"",$self->nt3->[$state],"\",(${tmp}*u,($yoffset+1.9)*u));\n";
		       print $mp "label.bot(\"",$codonAA[$state],"\",(${tmp}*u,$yoffset*u));\n"; }
    elsif ($N == 20) { print $mp "label.top(\"",$aa[$state],"\",(${tmp}*u,($yoffset+1.9)*u));\n"; }
  }

}


## Draw mutation rates between states ##
sub drawRates {
  my ($self, $mp, $dat, $showDiagonal) = @_;
  my $N = $self->N;
  my $penscale = $self->directives_penscale; my $rscale = $self->directives_rscale;

  if (!defined $showDiagonal) { $showDiagonal = 0; }

  #draw grid
  print $mp "% draw grid\n";
  if (!$self->directives_newgrid) { # a simple metapost loop to draw lines at every state
    print $mp "for i=1 upto $N:\n  draw (i*u,0)--(i*u,($N+1)*u) withpen pencircle scaled $penscale withcolor 0.6white;\n";
    print $mp "  draw (0,i*u)--(($N+1)*u,i*u) withpen pencircle scaled $penscale withcolor 0.7white;\nendfor;\n";
  }
  else { # metapost commands to draw lines between blocks of 4 codons with same 1st 2 nts, or between blocks of synonymous codons
    for (my $i = 0; $i < scalar(@{$self->grid}); $i++) {
      my $j = $N+1-$self->grid->[$i];
      print $mp "draw (",$self->grid->[$i],"*u,0)--(",$self->grid->[$i],"*u,($N+1)*u) withpen pencircle scaled $penscale withcolor 0.4white;\n";
      print $mp "draw (0,$j*u)--(($N+1)*u,$j*u) withpen pencircle scaled $penscale withcolor 0.4white;\n";
    }
  }

  #draw scale circle at top left if using a bubble scheme (i.e. not an experimental chess board scheme numbered >99)
  if ($self->directives_colsch < 100) {
    my $scalesize = 0.5;
    my $scaledvalue = 3 * $rscale * sqrt($scalesize);
    if ($self->directives_smooth) {
      print $mp "dullfullcircle((-2.4*u,-2.4*u),$scaledvalue*u,(0.8,0.8,0.8),$N);\n";
    }
    else {
      croak "prettyfullcircle stuff not implemented yet!\n";
#      print $mp "prettyfullcircle((-2.4*u,-2.4*u),$scaledvalue*u,$theta,$fraction,$smooth,(0.8,0.8,0.8),$N);\n";
    }
    print $mp "label(btex $scalesize etex scaled $rscale,(-2.4*u,($N+3.4)*u));\n";
  }

  my $count = 0; my (@x, @y, @z);
  #for each row of input file...
  while (<$dat>) {
    #pick out data
    my @no = split;
    if (scalar (@no) != 3) { croak ("Invalid input line '$_'."); }
    my ($fromState, $toState, $value) = ($no[0], $no[1], $no[2]);
    #if a valid matrix element...
    unless (($fromState < 1) || ($fromState > $N) || ($toState < 1) || ($toState > $N)) {

      # skip diagonal unless requested
      if ($fromState == $toState) { next unless $showDiagonal; }

      # skip stop codons
      if ($fromState > 61 || $toState > 61) { warn "Skipping transition $fromState -> $toState (stop codons disallowed)."; next; } # because we don't look at stop codons

      # have y-axis (row) be from state and x-axis (column) be to state
      my ($col, $row) = ($self->sorting->[$toState], $self->sorting->[$fromState]);
      if (!defined $col || !defined $row) { croak "Undefined column or row for input line '$_': from '$fromState' to '$toState'."; }

      if ($self->directives_colsch < 100) { #i.e. if a bubble scheme
	#output most of the circle command, except the color information
	my $scaledvalue = 3 * $rscale * sqrt($value);
	if ($self->directives_smooth) {
	  print $mp "dullfullcircle(($col*u,$row*u),";
	  printf $mp "%f*u,", $scaledvalue; # MetaPost seems to barf when $scaledvalue is in scientific notation; ugly printf hack to fix this
	}
	else {
	  croak "prettyfullcircle stuff not implemented yet!\n";
#	  print $mp "prettyfullcircle(($col*u,$row*u),";
#	  printf $mp "%f*u,", $scaledvalue;
#	  print $mp "$theta,$fraction,$smooth,";
	}
	#if needed, store values for re-drawing circles on top
	if ($self->directives_forceoutline) { $x[$count] = $col; $y[$count] = $row; $z[$count] = $scaledvalue; $count++; }
      }

      # this is where the coloring schemes start
      # each scheme computes the stats it needs,
      # and then completes the circle command output according to the programmed rules and colors

      if ($N == 4) {
	### scheme 0: color by transition/transversion
	if ($self->directives_colsch == 0) {
	  if (($pur[$fromState] && $pur[$toState]) || ((! $pur[$fromState]) && (! $pur[$toState]))) { print $mp "(0,1,0),$N);\n"; } # transition = green
	  else { print $mp "(0,0,0),$N);\n" } } # transversion
	### undefined scheme
	else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
      elsif ($N == 16) {
	### scheme 0: color by # nt differences
	if ($self->directives_colsch == 0) {
	  my $ndiffs = ($self->nt1->[$fromState] ne $self->nt1->[$toState])+($self->nt2->[$fromState] ne $self->nt2->[$toState]);
	  if    ($ndiffs == 1) { print $mp "(1,0,0),$N);\n"; } # 1 difference = red
	  elsif ($ndiffs == 2) { print $mp "(0,1,0),$N);\n"; } # 2 diffs = green
	  else                 { print $mp "(0,0,0),$N);\n"; } }  # else black
	### scheme 1: color by preservation of Watson-Crick pairing
	elsif ($self->directives_colsch == 1) {
	  if    ($wc[$fromState] && $wc[$toState])       { print $mp "(1,0,0),$N);\n"; } # Watson-Crick pairing preserved (WC -> WC) = red
	  elsif (($wc[$fromState]) && (! $wc[$toState])) { print $mp "(1,1,0),$N);\n"; } # WC -> non-WC = yellow
	  elsif ((! $wc[$fromState]) && ($wc[$toState])) { print $mp "(0,0,1),$N);\n"; } # non-WC -> WC = blue
	  else                                           { print $mp "(0,0,0),$N);\n" } } # non-WC -> non-WC = black
	### scheme 2: color by ts/tv
	elsif ($self->directives_colsch == 2) {
	  my $nts = (abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 2) + (abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 2); # number of transitions
	  my $ntv = (abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 1) + (abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 1); # number of A-C, C-G and G-T transversions...
	  $ntv += (abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 3)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 3); # ... plus number of A-T transversions
	  if    (($nts == 1) && ($ntv == 0)) { print $mp "(.941,0,0),$N);\n"; }          # dark red
	  elsif (($nts == 0) && ($ntv == 1)) { print $mp "(1,0.725,0.725),$N);\n"; }     # pale red
	  elsif (($nts == 2) && ($ntv == 0)) { print $mp "(0,0.6,0.047),$N);\n"; }       # dark green
	  elsif (($nts == 1) && ($ntv == 1)) { print $mp "(0.212,0.792,0.243),$N);\n"; } # mid green
	  elsif (($nts == 0) && ($ntv == 2)) { print $mp "(0.408,0.984,0.439),$N);\n"; } # pale green
	  else                               { print $mp "(0,0,0),$N);\n"; } }          # else black
	### undefined scheme
	else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
      elsif (($N == 61) || ($N == 64)) {
	### scheme 0: color by # nt differences
	if ($self->directives_colsch == 0) {
	  my $ndiffs = ($self->nt1->[$fromState] ne $self->nt1->[$toState])+($self->nt2->[$fromState] ne $self->nt2->[$toState])+($self->nt3->[$fromState] ne $self->nt3->[$toState]);
	  if    ($ndiffs == 1) { print $mp "(1,0,0),$N);\n"; } # 1 difference = red
	  elsif ($ndiffs == 2) { print $mp "(0,1,0),$N);\n"; } # 2 diffs = green
	  elsif ($ndiffs == 3) { print $mp "(0,0,1),$N);\n"; } # 3 diffs = blue
	  else                 { print $mp "(0,0,0),$N);\n"; } } # else black
	### scheme 1: color by synonymous/nonsynonymous
	elsif ($self->directives_colsch == 1) {
	  my $syn = ($codonAA[$fromState] eq $codonAA[$toState]);
	  if    ($syn == 1) { print $mp "(1,1,0),$N);\n"; } # synonymous = yellow
	  elsif ($syn == 0) { print $mp "(1,0,1),$N);\n"; } # nonsyn = purple
	  else              { print $mp "(0,0,0),$N);\n"; } } # else black
	### scheme 2: color by syn/nonsyn x 1 nt diff/>1 nt diff
	elsif ($self->directives_colsch == 2) {
	  my $ndiffs = ($self->nt1->[$fromState] ne $self->nt1->[$toState])+($self->nt2->[$fromState] ne $self->nt2->[$toState])+($self->nt3->[$fromState] ne $self->nt3->[$toState]);
	  my $syn = ($codonAA[$fromState] eq $codonAA[$toState]);
	  if    (($syn == 1) && ($ndiffs == 1)) { print $mp "(1,0,0),$N);\n"; } # syn & 1 diff = red
	  elsif (($syn == 1) && ($ndiffs>1))    { print $mp "(0,0,1),$N);\n"; } # syn & >1 diff = blue
	  elsif (($syn == 0) && ($ndiffs == 1)) { print $mp "(0,1,0),$N);\n"; } # nonsyn & 1 diff = green
	  elsif (($syn == 0) && ($ndiffs>1))    { print $mp "(1,1,0),$N);\n"; } # nonsyn & >1 diff = yellow
	  else                                  { print $mp "(0,0,0),$N);\n"; } } # else black
	### scheme 3: color by ts/tv
	elsif ($self->directives_colsch == 3) {
	  my $nts=(abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 2)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 2)+(abs($self->ntcode3->[$fromState]-$self->ntcode3->[$toState]) == 2); # number of transitions
	  my $ntv=(abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 1)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 1)+(abs($self->ntcode3->[$fromState]-$self->ntcode3->[$toState]) == 1); # number of A-C, C-G and G-T transversions...
	  $ntv+=(abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 3)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 3)+(abs($self->ntcode3->[$fromState]-$self->ntcode3->[$toState]) == 3); # ... plus number of A-T transversions
	  if    (($nts == 1) && ($ntv == 0)) { print $mp "(.941,0,0),$N);\n"; }          # dark red
	  elsif (($nts == 0) && ($ntv == 1)) { print $mp "(1,0.725,0.725),$N);\n"; }     # pale red
	  elsif (($nts == 2) && ($ntv == 0)) { print $mp "(0,0.6,0.047),$N);\n"; }       # dark green
	  elsif (($nts == 1) && ($ntv == 1)) { print $mp "(0.212,0.792,0.243),$N);\n"; } # mid green
	  elsif (($nts == 0) && ($ntv == 2)) { print $mp "(0.408,0.984,0.439),$N);\n"; } # pale green
	  elsif (($nts == 3) && ($ntv == 0)) { print $mp "(0,0,0.5),$N);\n"; }           # dark blue
	  elsif (($nts == 2) && ($ntv == 1)) { print $mp "(0.281,0.281,0.667),$N);\n"; } # mid-dark blue
	  elsif (($nts == 1) && ($ntv == 2)) { print $mp "(0.562,0.562,0.833),$N);\n"; } # mid-pale blue
	  elsif (($nts == 0) && ($ntv == 3)) { print $mp "(0.843,0.843,1),$N);\n"; }     # pale blue
	  else                               { print $mp "(0,0,0),$N);\n"; } }           # else black
	### scheme 4: CK's alternative coloring by ts/tv
	elsif ($self->directives_colsch == 4) {
	  my $nts=(abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 2)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 2)+(abs($self->ntcode3->[$fromState]-$self->ntcode3->[$toState]) == 2); # number of transitions
	  my $ntv=(abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 1)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 1)+(abs($self->ntcode3->[$fromState]-$self->ntcode3->[$toState]) == 1); # number of A-C, C-G and G-T transversions...
	  $ntv+=(abs($self->ntcode1->[$fromState]-$self->ntcode1->[$toState]) == 3)+(abs($self->ntcode2->[$fromState]-$self->ntcode2->[$toState]) == 3)+(abs($self->ntcode3->[$fromState]-$self->ntcode3->[$toState]) == 3); # ... plus number of A-T transversions
	  if    (($nts == 1) && ($ntv == 0)) { print $mp "(0.85,0.85,1),$N);\n"; } # light blue
	  elsif (($nts == 0) && ($ntv == 1)) { print $mp "(1,0.85,0.85),$N);\n"; } # light red
	  elsif (($nts == 2) && ($ntv == 0)) { print $mp "(0.5,0.5,1),$N);\n"; }   # medium blue
	  elsif (($nts == 1) && ($ntv == 1)) { print $mp "(1,0,1),$N);\n"; }       # medium violet
	  elsif (($nts == 0) && ($ntv == 2)) { print $mp "(1,0.5,0.5),$N);\n"; }   # medium red
	  elsif (($nts == 3) && ($ntv == 0)) { print $mp "(0,0,0.67),$N);\n"; }    # dark blue
	  elsif (($nts == 2) && ($ntv == 1)) { print $mp "(0.44,0,0.67),$N);\n"; } # dark blue-violet
	  elsif (($nts == 1) && ($ntv == 2)) { print $mp "(0.67,0,0.44),$N);\n"; } # dark red-violet
	  elsif (($nts == 0) && ($ntv == 3)) { print $mp "(0.67,0,0),$N);\n"; }    # dark red
	  else                               { print $mp "(0,0,0),$N);\n"; } }     # else black
	### scheme 5: color by # nt differences and which positions they are at
	elsif ($self->directives_colsch == 5) {
	  my $ndiffs_1=($self->nt1->[$fromState] ne $self->nt1->[$toState]);
	  my $ndiffs_2=($self->nt2->[$fromState] ne $self->nt2->[$toState]);
	  my $ndiffs_3=($self->nt3->[$fromState] ne $self->nt3->[$toState]);
	  my $ndiffs_tot=$ndiffs_1+$ndiffs_2+$ndiffs_3;
	  if    (($ndiffs_tot == 1) && ($ndiffs_2 == 1)) { print $mp "(1,0.8,0.8),$N);\n"; } # 1 diff,  @ pos 2   = pale red
	  elsif (($ndiffs_tot == 1) && ($ndiffs_1 == 1)) { print $mp "(1,0.4,0.4),$N);\n"; } # 1 diff,  @ pos 1   = med red
	  elsif (($ndiffs_tot == 1) && ($ndiffs_3 == 1)) { print $mp "(0.8,0,0),$N);\n"; }   # 1 diff,  @ pos 3   = dark red
	  elsif (($ndiffs_tot == 2) && ($ndiffs_3 == 0)) { print $mp "(0.8,1,0.8),$N);\n"; } # 2 diffs, @ pos 1+2 = pale green
	  elsif (($ndiffs_tot == 2) && ($ndiffs_1 == 0)) { print $mp "(0,0.9,0),$N);\n"; }   # 2 diffs, @ pos 2+3 = med green
	  elsif (($ndiffs_tot == 2) && ($ndiffs_2 == 0)) { print $mp "(0,0.6,0),$N);\n"; }   # 2 diffs, @ pos 1+3 = dark green
	  elsif  ($ndiffs_tot == 3)                      { print $mp "(0.8,0.8,1),$N);\n"; } # 3 diffs            = pale blue
	  else                                           { print $mp "(0,0,0),$N);\n"; } }  # else black
	### undefined scheme
	else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
      elsif ($N == 20) {
	my $fromAA = $aa[$fromState];
	my $toAA = $aa[$toState];
	### scheme 0
	if ($self->directives_colsch == 0) {
	  if ( grep($fromAA eq $_, @aaAromatic) && grep($toAA eq $_, @aaAromatic) ) {
	    print $mp "(1,0,0),$N);\n"; # red
	  }
	  elsif ( grep($fromAA eq $_, @aaHydrophobic) && grep($toAA eq $_, @aaHydrophobic) ) {
	    print $mp "(0,1,0),$N);\n"; # green
	  }
	  elsif ( grep($fromAA eq $_, @aaHydrophilic) && grep($toAA eq $_, @aaHydrophilic) ) {
	    print $mp "(0,0,1),$N);\n"; # blue
	  }
	  else {
	    print $mp "(0.5,0.5,0.5),$N);\n"; # grey
	  }	    
	}
	### undefined scheme
	else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
    }
  }

  #if needed, re-draw circles on top
  if (($self->directives_colsch < 100) && ($self->directives_forceoutline)) {
    for (my $i = 0; $i < $count; $i++) {
      print $mp "dullemptycircle(($x[$i]*u,$y[$i]*u),$z[$i]*u,$N);\n";
    }
  }

  # label axes
  if    ($N == 4)  { print $mp "label.lft(\"from\",(-u,($N+1)/2*u));\n";
		     print $mp "label.top(\"to\",(($N+1)/2*u,($N+2.0)*u));\n"; }
  elsif ($N == 16) { print $mp "label.lft(\"from\",(-2*u,($N+1)/2*u));\n";
		     print $mp "label.top(\"to\",(($N+1)/2*u,($N+3.4)*u));\n"; }
  elsif (($N == 61) || ($N == 64)) { print $mp "label.lft(\"from\",(-3*u,($N+1)/2*u));\n";
				     print $mp "label.top(\"to\",(($N+1)/2*u,($N+5)*u));\n"; }
  elsif ($N == 20) { print $mp "label.lft(\"from\",(-u,($N+1)/2*u));\n";
		     print $mp "label.top(\"to\",(($N+1)/2*u,($N+2.0)*u));\n"; }

  # label states
  for (my $toState = 1; $toState < $N+1; $toState++) {
    if ($toState > 61) { warn "Skipping codon '$toState' (stop codons disallowed)."; next; } # don't handle stop codons
    my $tmp = $N+1-$self->sorting->[$toState];
    if    ($N == 4)  { print $mp "label.lft(\"",$nuc[$toState],"\",(0,${tmp}*u));\n"; }
    elsif ($N == 16) { print $mp "label.lft(\"",$dinuc[$toState],"\",(0,${tmp}*u));\n"; }
    elsif (($N == 61) || ($N == 64)) {
		       print $mp "label.lft(\"",$codon[$toState],"\",(0,${tmp}*u));\n";
		       print $mp "label.rt(\"",$codonAA[$toState],"\",(($N+1)*u,${tmp}*u));\n"; }
    elsif ($N == 20) { print $mp "label.lft(\"",$aa[$toState],"\",(0,${tmp}*u));\n"; }
  }
  for (my $fromState = 1; $fromState < $N+1; $fromState++) {
    if ($fromState > 61) { warn "Skipping codon '$fromState' (stop codons disallowed)."; next; } # don't handle stop codons
    my $tmp = $self->sorting->[$fromState];
    if    ($N == 4)  { print $mp "label.top(\"",$nuc[$fromState],"\",(${tmp}*u,($N+1.0)*u));\n"; }
    elsif ($N == 16) { print $mp "label.top(\"",$self->nt1->[$fromState],"\",(${tmp}*u,($N+1.8)*u));\n";
		       print $mp "label.top(\"",$self->nt2->[$fromState],"\",(${tmp}*u,($N+1.0)*u));\n"; }
    elsif (($N == 61) || ($N == 64)) {
		       print $mp "label.top(\"",$self->nt1->[$fromState],"\",(${tmp}*u,($N+2.7)*u));\n";
		       print $mp "label.top(\"",$self->nt2->[$fromState],"\",(${tmp}*u,($N+1.9)*u));\n";
		       print $mp "label.top(\"",$self->nt3->[$fromState],"\",(${tmp}*u,($N+1.1)*u));\n";
		       print $mp "label.bot(\"",$codonAA[$fromState],"\",(${tmp}*u,0));\n"; }
    elsif ($N == 20) { print $mp "label.top(\"",$aa[$fromState],"\",(${tmp}*u,($N+1.0)*u));\n"; }
  }

}



###########################
## Graphics-related subs ##
###########################

sub startFig {
  my ($self, $mp) = @_;

  print $mp "%%%%%%%%%%\n
beginfig(1)\n
defaultscale := 0.5;
pickup pencircle scaled 0.5pt;\n
u := 6pt;\n
%%%%%%%%%%\n\n";
}

sub endFig {
  my ($self, $mp) = @_;

  print $mp "\n\nendfig;\nend;\n";
}

sub drawInitBox {
  my ($self, $mp) = @_;
  my $N = $self->N;
  my ($xoffset, $yoffset) = ($self->directives_xoffset,$self->directives_yoffset);

  print $mp "\n\n%%%%%%%%%%\n\n";
  print $mp "% draw outline box for initial probability distribution\n";
  print $mp "draw ((-1+$xoffset)*u,$yoffset*u)--(($N+2+$xoffset)*u,$yoffset*u);\n"; # top
  print $mp "draw ((0+$xoffset)*u,($yoffset-1)*u)--((0+$xoffset)*u,($yoffset+3)*u);\n"; # left
  print $mp "draw ((-1+$xoffset)*u,($yoffset+2)*u)--(($N+2+$xoffset)*u,($yoffset+2)*u);\n"; # bottom
  print $mp "draw (($N+1+$xoffset)*u,($yoffset-1)*u)--(($N+1+$xoffset)*u,($yoffset+3)*u);\n\n"; # right
  print $mp "% tickmarks\n";
  print $mp "draw ((-1/4+$xoffset)*u,($yoffset+1)*u)--((1/4+$xoffset)*u,($yoffset+1)*u);\n"; # left
  print $mp "draw (($N+1+$xoffset)*u-u/4,($yoffset+1)*u)--(($N+1+$xoffset)*u+u/4,($yoffset+1)*u);\n"; # right
  print $mp "for i=1 upto $N:\n";
  print $mp "  draw ((i+$xoffset)*u,($yoffset+1/4)*u)--((i+$xoffset)*u,($yoffset-1/4)*u);\n"; # top
  print $mp "  draw ((i+$xoffset)*u,($yoffset+2+1/4)*u)--((i+$xoffset)*u,($yoffset+2-1/4)*u);\n"; # bottom
  print $mp "endfor;\n\n";
}

sub drawRatesBox {
  my ($self, $mp) = @_;
  my $N = $self->N;

  print $mp "\n\n%%%%%%%%%%\n\n";
  print $mp "% draw main outline box for mutation rates\n";
  print $mp "draw (-u,0)--(($N+2)*u,0);\ndraw (0,-u)--(0,($N+2)*u);\ndraw (-u,($N+1)*u)--(($N+2)*u,($N+1)*u);\ndraw (($N+1)*u,-u)--(($N+1)*u,($N+2)*u);\n\n";
  print $mp "% tickmarks\n";
  print $mp "for i=1 upto $N:\n";
  print $mp "  draw (i*u,($N+1)*u+u/4)--(i*u,($N+1)*u-u/4);\n";
  print $mp "  draw (-u/4,($N+1)*u-i*u)--(u/4,($N+1)*u-i*u);\n";
  print $mp "  draw (($N+1)*u-u/4,($N+1)*u-i*u)--(($N+1)*u+u/4,($N+1)*u-i*u);\n";
  print $mp "  draw (i*u,u/4)--(i*u,-u/4);\n";
  print $mp "endfor;\n\n";
}

sub printTex {
  my ($self, $tex, $files, $legends, $captionFile) = @_; # $legends is an array reference (legend for each file in @files)
  my $figscale = $self->directives_figscale;

  print $tex "\\documentclass[english]{article}
\\usepackage[T1]{fontenc}
\\usepackage[latin1]{inputenc}
\\usepackage{graphicx}
\\usepackage{caption}
\\usepackage{color}

\\makeatletter
\\usepackage{babel}
\\makeatother

\\pagestyle{empty}

\\begin{document}\n";

  foreach my $file (@$files) {
    my $legend = shift @$legends;
    print $tex "\\begin{figure}
\\centering
\\includegraphics[scale=$figscale]{./$file}\n";
    if ($captionFile) {
      print $tex "\\input{$captionFile}\n";
    }
    else {
      # stupid bugfix for empty legend... because apparently the below doesn't work on empty legends if
      # you have an incomplete or outdated TeX (?) installation, so we should bypass a potential error in
      # cases where a legend in the caption is simply not wanted
      unless ($legend eq '') {
	print $tex "\\captionsetup{labelformat=empty, labelsep=none, aboveskip=0cm, textfont=footnotesize, justification=centering}
\\caption{$legend}\n";
      }
    }
    print $tex "\\end{figure}\n"; 
  }

  print $tex "\\end{document}\n";
}

# stub of a function to automatically choose a default figure legend based on the alphabet size and color scheme
sub defaultLegend {
  my ($self) = @_;
  my $N = $self->N;
  my $legend = "No default legend defined!";

  if ($N == 4) {
    ### scheme 0: color by transition/transversion
    if ($self->directives_colsch == 0) { }
    ### undefined scheme
    else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
  elsif ($N == 16) {
    ### scheme 0: color by # nt differences
    if ($self->directives_colsch == 0) { }
    ### scheme 1: color by preservation of Watson-Crick pairing
    elsif ($self->directives_colsch == 1) { $legend = "\$\\qquad\\qquad\\,\$(Red, Yellow): WC \$\\longrightarrow\$ (WC, non-WC) \\newline (Blue, Black): non-WC \$\\longrightarrow\$ (WC, non-WC)"; }
    ### scheme 2: color by ts/tv
    elsif ($self->directives_colsch == 2) { }
    ### undefined scheme
    else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
  elsif (($N == 61) || ($N == 64)) {
    ### scheme 0: color by # nt differences
    if ($self->directives_colsch == 0) { }
    ### scheme 1: color by synonymous/nonsynonymous
    elsif ($self->directives_colsch == 1) { }
    ### scheme 2: color by syn/nonsyn x 1 nt diff/>1 nt diff
    elsif ($self->directives_colsch == 2) { }
    ### scheme 3: color by ts/tv
    elsif ($self->directives_colsch == 3) { }
    ### scheme 4: CK's alternative coloring by ts/tv
    elsif ($self->directives_colsch == 4) { }
    ### scheme 5: color by # nt differences and which positions they are at
    elsif ($self->directives_colsch == 5) { }
    ### undefined scheme
    else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }
  elsif ($N == 20) {
    ### scheme 0: color all same
    if ($self->directives_colsch == 0) { 
      $legend = "\\textcolor{red}{Red}: Aromatic \$\\longrightarrow\$ Aromatic \\newline " . 
	  "\\textcolor{green}{Green}: Hydrophobic \$\\longrightarrow\$ Hydrophobic \\newline " .
	  "\\textcolor{blue}{Blue}: Hydrophillic \$\\longrightarrow\$ Hydrophillic \\newline " .
          "Gray: Any of the above \$\\longrightarrow\$ Another group \\newline "; 
   }
    ### undefined scheme
    else { croak "\n\nWARNING: Color scheme ",$self->directives_colsch," not defined for N = $N\n\n\n"; } }

  return $legend;
}


sub printDefs {
  my ($self, $mp) = @_;
  my $penscale = $self->directives_penscale;

  print $mp "input rboxes;

def prettyfullcircle(expr c,r,t,l,n,B,N) =
% this could really use some documentation!
% routine abandoned because it caused memory problems  
%
% flip y-coordinate so that we're labelling from top left
  x:=xpart c; y:=(N+1)*u-(ypart c);
%
  for i=0 upto n:
%    show ((n-i)*2r/n),((x,y)+(i*l*r/n)*dir(t)),((n-0.9*i)*B/n+(0.9*i)*white/n);
    fill fullcircle scaled ((n-i)*2r/n) shifted ((x,y)+(i*l*r/n)*dir(t)) withcolor ((n-0.9*i)*B/n+(0.9*i)*white/n);
  endfor;
  draw fullcircle scaled 2r shifted (x,y) withpen pencircle scaled $penscale;
enddef;

def dullsquare(expr c,B) =
% draw a coloured square of 'unit' size centred on c with colour B
%
% flip y-coordinate so that we're labelling from top left
  x:=xpart c; y:=17u-(ypart c);
%
  path p;
  p=(x-u/2,y-u/2)--(x-u/2,y+u/2)--(x+u/2,y+u/2)--(x+u/2,y-u/2)--cycle;
%  draw p withpen pencircle scaled $penscale;
  fill p withcolor B;
%
enddef;

def dullfullcircle(expr c,r,B,N) =
% draw a coloured disk at c, radius r, colour B, with black outline
% called 'dull' because it's not as pretty as prettyfullcircle above!
% N is dimension of state space + 1
% flip y-coordinate so that we're labelling from top left
  x:=xpart c; y:=(N+1)*u-(ypart c);
%
% draw disk and outline, or just a grey dot if radius is 0
  if r>0 :
    fill fullcircle scaled 2r shifted (x,y) withcolor B;
    draw fullcircle scaled 2r shifted (x,y) withpen pencircle scaled $penscale;
  else :
    draw fullcircle scaled 0 shifted (x,y) withpen pencircle scaled $penscale withcolor 0.7white;
  fi;
%
enddef;

def dullemptycircle(expr c,r,N) =
% draw a circle at c, radius r
%
% flip y-coordinate so that we're labelling from top left
  x:=xpart c; y:=(N+1)*u-(ypart c);
%
% draw circle, or just a grey dot if radius is 0
  if r>0 :
    draw fullcircle scaled 2r shifted (x,y) withpen pencircle scaled $penscale;
  else :
    draw fullcircle scaled 0 shifted (x,y) withpen pencircle scaled $penscale withcolor 0.7white;
  fi;
%
enddef;\n\n";
}



1
