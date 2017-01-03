#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

=head1 NAME

phylogramTool.pl - Various methods for operating on phylograms.

=cut

use strict;
use Getopt::Long;
use Phylogram;

# Read command line arguments

my ($help, $expChainRate, $scaleFactor);

GetOptions( "h"    => \$help,
            "ecr"  => \$expChainRate,
            "s=s"  => \$scaleFactor);
	     
my $pgFile = shift;

if( $help or not $pgFile) {
    help();
    exit(1);
}

# Main
my $pg=PhyloGram->from_file($pgFile);
if ($expChainRate)
{
  calcExpChainRate($pg);
}
elsif ($scaleFactor)
{
  scaleRates($pg, $scaleFactor);
}

# Subroutines

=head1 METHOD
calcExpChainRate - calculate expected rate of substitution for each chain in phylogrammar
Input: 
  pg - phylogrammar
Return:
  None
=cut
sub calcExpChainRate
{
  my $pg = shift();
  foreach my $chain ($pg->all_chains)
  {
    print $chain->terminal->to_string,"\n";
    my $mutateHash=$chain->mutate_hash;
    my @states=map(join("",@{$_->state->value}), $chain->find_all("initial"));
    print"States: @states\n";
    my $total=0;
    foreach my $i (@states)
    {
      foreach my $j (@states)
      {
	if ($i ne $j)
	{
	  my $m=$chain->mutate($i,$j,$mutateHash);
	  if(defined $m)
	  {
	    $total+=$chain->initial($i)->prob->value * $m->rate->value
	  }
	}
      }
    }
    print"Expected mutation rate: $total\n","-"x80,"\n"
  }
}

=head1 METHOD
scaleRates - scale the phylogrammar rates by a constant and output.
Input: 
  pg - phylogrammar
  mul - scaling factor
Return:
  None
=cut
sub scaleRates
{
  my $pg = shift();
  my $mul=shift();
  foreach my $chain ($pg->all_chains)
  {
    my $mutateHash=$chain->mutate_hash;
    my @states=map(join("",@{$_->state->value}),$chain->find_all("initial"));
    foreach my $i (@states)
    {
      foreach my $j (@states)
      {
	if ($i ne $j) 
	{
	  my $m=$chain->mutate($i,$j,$mutateHash);
	  if(defined $m)
	  {
	    $m->rate->value ($m->rate->value * $mul)
	  }
	}
      }
    }
  }
  print $pg->to_string;
}

sub help {
    print STDERR <<EOF;

$0: 
Various methods for operating on phylograms.

Usage: $0 [options] phylogram
    Options
        -h             : show this help
        -ecr           : Calculate expected substitution rate of every chain in grammar.
        -s num         : Scale phylogram rates by a constant.  
                         Use xgram -t new.eg -g old.eg < /dev/null to reformat.
EOF
}
