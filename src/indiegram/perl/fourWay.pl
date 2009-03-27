#!/usr/bin/perl -w

use ComposedTreeTransducer::FourWayComposedTT;

use strict;
use Carp;

# NB: Delete states are special cases of 
# match states and should be implemented as such.

# Bifurcations also have emit and absorb profiles.
# Bifurcations are deterministic (probability 1).
# Child states must be of type Start.

# The (unique) "overall" start states of the singlet and branch 
# machines are taken to be the first states in the transition matrix

my $usage = "\nUsage: $0 <.singlet.tt file> <.branch.tt file>
           [-h, --help] display this message
           [-d, --debug] display debugging info
           [-r, --redundant] look for redundant states in the reduced TM (costly!)
           [--nosuppress] (default is to suppress) don't suppress all-End children of bifurcation states
           [--noprune] (default is to prune) don't remove Bifurc states with only 1 non-End child
           [--graph] display state graph
           [--tm] display transition matrix
           [--rgraph] display reduced state graph
           [--rtm] display reduced transition matrix
           [--frgraph] display fully reduced state graph
           [--frtm] display fully reduced transition matrix
           [-a, --all] display all the info I've got
           [--tkfst] use the TKFST modelj (uses ../tt/tkfst.{singlet,branch}.tt)

           For information on the .tt file format, see http://biowiki.org/TreeTransducerComposition.\n";

my ($singletfile, $branchfile);
my $directives = {}; # hash to store command-line flags

while (@ARGV) {
  my $arg = shift @ARGV;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif (($arg eq "-d") || ($arg eq "--debug")) { $directives->{'debug'} = 1; }
    elsif (($arg eq "-r") || ($arg eq "--redundant")) { $directives->{'redundant'} = 1; }
    elsif (($arg eq "--nosuppress")) { $directives->{'nosuppress'} = 1; }
    elsif (($arg eq "--noprune")) { $directives->{'noprune'} = 1; }
    elsif ($arg eq "--graph") { $directives->{'graph'} = 1; }
    elsif (($arg eq "--tm")) { $directives->{'tm'} = 1; }
    elsif (($arg eq "--rgraph")) { $directives->{'rgraph'} = 1; }
    elsif (($arg eq "--rtm")) { $directives->{'rtm'} = 1; }
    elsif (($arg eq "--frgraph")) { $directives->{'frgraph'} = 1; }
    elsif (($arg eq "--frtm")) { $directives->{'frtm'} = 1; }
    elsif (($arg eq "-a" ) || ($arg eq "--all")) { $directives->{'all'} = 1; }
    elsif ($arg eq "--tkfst") { $singletfile = "../tt/tkfst.singlet.tt"; $branchfile = "../tt/tkfst.branch.tt"; }
    else { croak $usage; }
  }
  elsif ($arg =~ /\S+\.singlet\.tt$/) { $singletfile = $arg; }
  elsif ($arg =~ /\S+\.branch\.tt$/) { $branchfile = $arg; }
  else { croak $usage; }
}

unless (defined $singletfile && defined $branchfile) { print $usage; exit; }
unless (-e ($singletfile)) { croak ("'$singletfile' not found.\n"); }
unless (-e ($branchfile)) { croak ("'$branchfile' not found.\n"); }

##########
# Get TM #
##########

my $composedTT = ComposedTreeTransducer::FourWayComposedTT->new ($singletfile, $branchfile, $directives);

$composedTT->show();
