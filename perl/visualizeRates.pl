#!/usr/bin/env perl -w

use Carp;

use PhyloGram;
use ChainDat;
use Bubbles;

my $usage = "\nUsage: $0 <list of files to process ('.dat' or xgram '.eg' grammar files)> [options]

           Name of output file = <name of input file>_<chain #>, with extention changed to '.ps' or '.eps'.

           For '.dat' files, please merge both rates and initial probabilities into the same file - the program will figure out which is which.
           Nucleotides, dinucleotides, amino acids and codons (stop codons removed) must be numbered (1-based indexing) alphabetically in '.dat.' files.

           [-h, --help] display this message
           [-c, --colorscheme] (nuc)   0 = ts/tv
                               (dinuc) 0 = # nt diffs (default); 1 = preservation of canonical pairing; 2 = ts/tv
                               (codon) 0 = # nt diffs (default); 1 = synonymous/nonsynonymous;
                                       2 = color by syn/nonsyn x 1 nt diff/>1 nt diff
                                       3 = ts/tv; 4 = alternate ts/tv; 5 = # nt differences and which positions they are at
                               (a.a.)  0 = all same (default)
           [-cap, --caption] caption tex file to be included in figure
           [-ch, --chain] chain number
           [-scale, --figscale] scale factor for the entire bubbleplot (default varies with alphabet)
           [-rscale, --ratescalefactor] rate bubble size scale factor (default is 1)
           [-pscale, --probscalefactor] initial probability distribution scale factor (default varies with alphabet)
           [-r, --resort] (codon) sort codons by: 0 = alphabetical (default), 1 = genetic ordering, 2 = a.a. identity, 3 = biochemical similarity, 4 = sorting due to Urbina and Tang (JME 2006) and Higgs
              (((All of the above flags take numerical arguments.  Below flags take no argument.)))
           [--newgrid] (codon) grid by synonymous blocks (default is lines at every codon)
	   [-oc] display observed-chain-counts instead of rate matrix for selected chain (experimental)
           [--forceoutline] redraw circles on top
           [--ratesonly] display only mutation rates (default is to print initial probability distribution as well)
           [--nolegend] no figure legend
           [--toeps] output file in Encapsulated Postscript format
           [--rnamode] convert all T bases to U's for output
           [--gupairs] treat GU base pairs as canonical base pairs
\nHandles 1st, 2nd and 3rd order Markov chains over nucleotides (RNA or DNA alphabet), corresponding respectively
to single nuc, dinuc and codon evolution, or 1st order chains over amino acids.
See http://biowiki.org/BubblePlots for more information and some tips.\n\n";

my $directives; # hash to store command-line flags
my ($inputChainIndex, $captionFile);
my (@infiles);
my ($nolegend, $toeps);

while (@ARGV) {
  my $arg = shift @ARGV;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit }
    elsif (($arg eq "-c") || ($arg eq "--colorscheme")) { $directives->{'colsch'} = shift @ARGV; }
    elsif (($arg eq "-cap") || ($arg eq "--caption")) { $captionFile = shift @ARGV; }
    elsif (($arg eq "-ch") || ($arg eq "--chain")) { $inputChainIndex = shift @ARGV; }
    elsif (($arg eq "-scale") || ($arg eq "--figscale")) { $directives->{'figscale'} = shift @ARGV; }
    elsif (($arg eq "-rscale") || ($arg eq "--ratescalefactor")) { $directives->{'rscale'} = shift @ARGV; }
    elsif (($arg eq "-pscale") || ($arg eq "--probscalefactor")) { $directives->{'pscale'} = shift @ARGV; }
    elsif (($arg eq "-r") || ($arg eq "--resort")) { $directives->{'resort'} = shift @ARGV; }
    elsif ($arg eq "--newgrid") { $directives->{'newgrid'} = 1; }
    elsif ($arg eq "--forceoutline") { $directives->{'forceoutline'} = 1; }
    elsif ($arg eq "--ratesonly") { $directives->{'ratesonly'} = 1; }
    elsif ($arg eq "--nolegend") { $nolegend = 1; }
    elsif ($arg eq "--toeps") { $toeps = 1; }
    elsif ($arg eq "--rnamode") { $directives->{'rnamode'} = 1; }
    elsif ($arg eq "--gupairs") { $directives->{'gupairs'} = 1; } 
    elsif ($arg eq "-oc") { $directives->{'oc'} = 1; }
    else { croak $usage; }
  }
#  elsif ($arg =~ /\.ps$/) { $outfile = $arg; } # deprecated functionality
#  elsif ($arg =~ /\.eps$/) { $outfile = $arg; $toeps = 1; }
  elsif ($arg =~ /\.eg$/) { push(@infiles, $arg); }
  elsif ($arg =~ /\.dat$/) { push(@infiles, $arg); }
  else { croak $usage; }
}

unless (@infiles > 0) { print $usage; exit; }

if (defined($captionFile) && ! -e ($captionFile)) {
  croak("Caption file [$captionFile] not found.\n");
}

############ Start actual script  ############

# Deal with command-line flags
unless (defined $directives->{'colsch'}) { $directives->{'colsch'} = 0; }
unless (defined $directives->{'figscale'}) { $directives->{'figscale'} = 1; }

$directives->{'smooth'} = 1; # force usage of dullfullcircle

unless (defined $directives->{'resort'}) { $directives->{'resort'} = 0; }

# Iterate over all the input files
foreach my $infile (@infiles) {
  my ($dat_rates_file, $dat_init_file, $N, $outfile);

  # Shall we load from a .dat file or an xgram file?

  if ($infile =~ /(.*)\.dat$/) {  # loading rates from .dat file
    $dat_rates_file = $1 . '_tmprates.dat';
    $dat_init_file = $1 . '_tmpinit.dat';

    # name output file
    if ($toeps) { $outfile = $1 . '.eps'; }
    else        { $outfile = $1 . '.ps'; }
    $outfile =~ s/^.*\///;  # make the output file local

    # Sanitize input .dat file, copying & converting to a clean form; count the number of states

    my %states;  # for storing all unique state numbers that we encounter in the input file
    my @rates; # for automatic sizing

    open (INDAT, "<$infile");
    open (OUTDAT_RATES, ">$dat_rates_file");
    open (OUTDAT_INIT, ">$dat_init_file");

    my $linenum = 1;
    while (my $line = <INDAT>) {
      unless ($line =~ /^\s+$/) {  # filter out blank lines
	my @row = split (' ', $line);

	if (@row == 3) {  # row contains 3 columns; assume it is a rate matrix row ('from' 'to' 'rate')
 	  $states{$row[0]} = 1;
	  $states{$row[1]} = 1;
	  push @rates, $row[2];
	  print OUTDAT_RATES $line;
	}
	elsif (@row == 2) {  # row contains 2 columns; assume it is an initial probability row
	  $states{$row[0]} = 1;
	  print OUTDAT_INIT $line;
	}
	else {
	  warn "ERROR: invalid entry on line $linenum of input file $infile - skipping...\n";
	}
      }
      $linenum++;

    }

    # automatic default sizing
    unless (defined $directives->{'rscale'}) {
      @rates = sort {$a < $b} @rates;
      #auto rate formulas totally arbitrary, change at will
      if ($rates[0] > 1 ) {$directives->{'rscale'} = 1/$rates[0];}
      else {$directives->{'rscale'} = .1 + 1-$rates[0];}
    }

    close (INDAT);
    close (OUTDAT_RATES);
    close (OUTDAT_INIT);

    $N = keys %states;  # now, the number of unique states should be the number of keys in the hash
  }
  elsif ($infile =~ /(.*)\.eg$/) {  # loading rates from an xrate grammar file
    $dat_rates_file = $1 . '_tmprates.dat';
    $dat_init_file = $1 . '_tmpinit.dat';

    my $grammar = PhyloGram->from_file($infile);
   
    # Select a chain (PhyloGram method all_chains returns a list of objects of type PhyloGram::Chain)
    unless (defined $inputChainIndex) { $inputChainIndex = selectChain ($grammar->all_chains); }
    my $chain = ($grammar->all_chains())[$inputChainIndex];

    # name output file
    if ($toeps) { $outfile = "$1\_$inputChainIndex" . '.eps'; }
    else        { $outfile = "$1\_$inputChainIndex" . '.ps'; }
    $outfile =~ s/^.*\///;  # make the output file local

    # read in parameters from grammar file
    my $params;
    my $h = $grammar->param_hash;
    while(($n,$s) = each %$h) { $params->{$n} = $s->value; }

    if (defined $directives->{'oc'}) {
      $oc = 'observed-chain-counts';
      @ochains = $grammar->grammar->$oc->values;
      @ochain = $ochains[$inputChainIndex]->values;
      my $high_init = 0;
      my $high_mutate = 0;
      foreach $piece (@ochain){
	if($piece->tag eq 'initial'){
	  $high_init += $piece->count->value;
	} elsif($piece->tag eq 'mutate'){
	  if ($piece->count->value > $high_mutate){
	    $high_mutate = $piece->count->value;
	  }
	}
      }
      
      foreach $piece (@ochain){
	if($piece->tag eq 'initial'){
	  $state = $piece->state->value;
	  $prob = $piece->count->value / $high_init;
	  $chain->initial($state, $prob);
	}elsif($piece->tag eq 'mutate'){
	  $from = $piece->from->value;
	  $to = $piece->to->value;
	  $rate = $piece->count->value / $high_mutate;
	  $chain->mutate($from, $to, $rate);
	}
      }
      
    }
       
    $chaindat = ChainDat->new($params, $chain);
    $N = $chaindat->N;

    # dump rates and initial probabilities to their respective temporary .dat files
    open DAT, ">$dat_rates_file" or croak "Couldn't create '$dat_rates_file': $!";
    $chaindat->ratesToDat(\*DAT);
    close DAT;

    open DAT, ">$dat_init_file" or croak "Couldn't create '$dat_init_file': $!";
    $chaindat->initToDat(\*DAT);
    close DAT;

    #automatic default sizing, w00t! - LEB 12/17/2007
    unless (defined $directives->{'rscale'}) {
      my @mutates = $chain->find_all("mutate");
      my @rates;
      foreach (@mutates) {
	@values = $_->rate->values;
	push(@rates, $chaindat->parseParams(\@values));
      }
      @rates = sort {$a < $b} @rates;
      #auto rate formulas totally arbitrary, change at will
      if($rates[0] > 1 ){$directives->{'rscale'} = 1/$rates[0];}
      else {$directives->{'rscale'} = .1 + 1-$rates[0];}
    }
  }
  else {
    # just in case...
    die ("ERROR: unsupported file format: '$infile' (but, this should not have occured, it's a bug).\n");
  }
	
  ## The rest of the code should be the same regardless of whether rates came from .dat or xrate file ## 
  # Do alphabet-size dependent calculations:
  #  resize bubbles in initial probability distribution plot
  unless (defined $directives->{'pscale'}) {
    if    ($N == 4)                  { $directives->{'pscale'} = 0.3; }
    elsif ($N == 16)                 { $directives->{'pscale'} = 0.8; }
    elsif (($N == 61) || ($N == 64)) { $directives->{'pscale'} = 1.5; }
    elsif ($N == 20)                 { $directives->{'pscale'} = 1; }
  }
  #  resize figures in latex file for aesthetics
  #  make sure grid lines and circle borders scale properly with figures
  $directives->{'penscale'} = 0.1; # original value in Nick's program
  if    ($N == 4)                  { $directives->{'figscale'} *= 3; $directives->{'penscale'} *= 6; }
  elsif ($N == 16)                 { $directives->{'figscale'} *= 2; $directives->{'penscale'} *= 4; }
  elsif (($N == 61) || ($N == 64)) { $directives->{'figscale'} *= 1; $directives->{'penscale'} *= 1; }
  elsif ($N == 20)                 { $directives->{'figscale'} *= 1.8; $directives->{'penscale'} *= 3.5; }

  # Set offset of plot of initial probability distribution relative to origin (used to get it properly aligned with mutation rate plot)
  if (($N == 61) || ($N == 64))    { ($directives->{'xoffset'},$directives->{'yoffset'}) = (0,-7); }
  else                             { ($directives->{'xoffset'},$directives->{'yoffset'}) = (0,-5); }

  ## Start figure in Metapost file ##
  open MP, ">bubbles.mp" or croak "Couldn't open 'bubbles.mp': $!";
  my $figure = Bubbles->new($directives, $N);
  $figure->printDefs(\*MP);
  $figure->startFig(\*MP);

  system "sort -rk3 $dat_rates_file > tmp";
  system "mv tmp $dat_rates_file";
  open DAT, "<$dat_rates_file" or croak "Couldn't open '$dat_rates_file': $!";
  $figure->drawRatesBox(\*MP);
  $figure->drawRates(\*MP, \*DAT);
  close DAT;
  system "rm $dat_rates_file";

  unless ($figure->directives_ratesonly) { # unless only printing rates
    #put the bubbles into 'biggest first' order
    system "sort -rk3 $dat_init_file > tmp";
    system "mv tmp $dat_init_file";
    open DAT, "<$dat_init_file" or croak "Couldn't open '$dat_init_file': $!";
    $figure->drawInitBox(\*MP);
    $figure->drawInitDist(\*MP, \*DAT);
    close DAT;
    system "rm $dat_init_file";
  }

  $figure->endFig(\*MP);
  close MP;
  system "mpost bubbles.mp";

  # Set figure legend
  #  standard latex math directives are ok
  #  the method defaultLegend() is a stub which only knows how to handle an alphabet of size 16 with colorscheme 1; consult it for an
  #   example legend
  my $legend = "";
  unless ($nolegend) { $legend = $figure->defaultLegend(); }

  ## Make final figure
  open TEX, ">bubbles.tex" or croak "Couldn't open 'bubbles.tex': $!\n";
  $figure->printTex(\*TEX, ["bubbles.1"], [$legend], $captionFile);
  close TEX;
  system "latex bubbles";
  system "dvips -o bubbles.ps bubbles";
  if ($toeps) {
    print "\nDoing ps2epsi...\n";
    system "ps2epsi bubbles.ps bubbles.epsi";
    print "Doing eps2eps...\n";
    system "eps2eps bubbles.epsi bubbles.eps";
    system "mv bubbles.eps $outfile";
  } else {
    system "mv bubbles.ps $outfile";
  }
  print "\nCreated '$outfile'\n\n";

  # clean up
  system "rm bubbles.mp bubbles.tex bubbles.1 bubbles.dvi bubbles.aux bubbles.log bubbles.mpx";
}


##########
## Subs ##
##########

## Select chain with user input.  Returns index of chain.
sub selectChain {
  my @chains = @_;
  my $chain;

  print "\n",scalar(@chains)," chains available:\n";
  for (my $k = 0; $k < scalar(@chains); $k++) {
    $chain = $chains[$k];
    print "\n($k).\n";
    print $chain->terminal->to_string, "\n"; # print out the terminals, using DartSexpr's AUTOLOADed methods to find the child node of "chain" named "terminal" ("to_string" is a built-in DartSexpr method)
    my @states = map(join("", @{$_->state->value}), $chain->find_all("initial"));
    print "states: @states\n";
  }

  if (@chains < 2) { print "\n(auto-selecting chain (0))\n\n"; return 0; }

  print "\nSelect a chain (0-",scalar(@chains)-1,"): ";
  $k = <STDIN>; chomp($k);
  if (($k >= 0) && ($k < scalar(@chains))) {
    return $k;
  } else { croak "Invalid selection.  Only chains 0-",scalar(@chains)-1," available!\n"; }

}

# MISC old code

#  if (($value =~ /^\d*.?\d*$/g) || ($value =~ /^\d*.?\d*e-\d+$/g)) { print $dat "$value\n"; } # if decimal or scientific notation number
#  else {
#  $value =~ /\D*(\d*).(\d*)\D*/g; # strip off parametric info
#  my $newvalue = $1 . "." . $2;
#  if (!($value eq $newvalue)) { print "Parametric info removed: $value -> $newvalue\n"; $parametricFlag = 1; }
#  print $dat "$newvalue\n"; }
