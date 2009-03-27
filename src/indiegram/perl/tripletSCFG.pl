#!/usr/bin/perl -w

use ComposedTreeTransducer::TripletSCFG;

use strict;
use Carp;

my $usage = "\nUsage: $0 <.singlet.tt file> <.branch.tt file>
           [-h, --help] display this message
           [-d, --debug] display debugging info
           [-r, --redundant] look for redundant states in the reduced TM (costly!)
           [--tkfst] build the TKFST_Triplet_SCFG (uses ../tt/tkfst.{singlet,branch}.tt)
           [--write] display the formatted TripletSCFG

           For information on the .tt file format, see http://biowiki.org/IndiegramSoftware.
           The purpose of this script is to coerce the output of a four-way tree transducer composition
           (see TreeTransducer::FourWayComposedTT) into the format needed for my Triplet_SCFG code.\n";

my ($singletfile, $branchfile);
my $directives = {}; # hash to store command-line flags

# set defaults
my $debug = 0;
my $nosuppress = 0;

while (@ARGV) {
  my $arg = shift @ARGV;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif (($arg eq "-r") || ($arg eq "--redundant")) { $directives->{'redundant'} = 1; }
    elsif ($arg eq "--tkfst") { $singletfile = "../tt/tkfst.singlet.tt"; $branchfile = "../tt/tkfst.branch.tt"; }
    elsif (($arg eq "--write")) { $directives->{'write'} = 1; }
    else { croak $usage; }
  }
  elsif ($arg =~ /\S+\.singlet\.tt$/) { $singletfile = $arg; }
  elsif ($arg =~ /\S+\.branch\.tt$/) { $branchfile = $arg; }
  else { croak $usage; }
}

unless (defined $singletfile && defined $branchfile) { print $usage; exit; }
unless (-e ($singletfile)) { croak ("'$singletfile' not found.\n"); }
unless (-e ($branchfile)) { croak ("'$branchfile' not found.\n"); }

#############
# Do Stuff! #
#############

my $triplet = ComposedTreeTransducer::TripletSCFG->new ($singletfile, $branchfile, $directives);

$triplet->show();
