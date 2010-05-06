#!/usr/bin/env perl -w

=head1 NAME

stkCompareAnnot.pl - Compares the column annotations of two stockholm alignments.

=cut

use strict;
use Getopt::Long;

use Roc;
use Stockholm::Database;

#my $mapDssp3 = 0;
#&GetOptions( "mapDssp3" => \$mapDssp3);
if (@ARGV != 3) { help(); exit(1)}
my $inputFeature = shift;
my $dbfile1 = shift;
my $dbfile2 = shift;

my $dbArRef1 = Stockholm::Database->from_file ($dbfile1);
my $dbArRef2 = Stockholm::Database->from_file ($dbfile2);
my $db1 = @$dbArRef1[0];
my $db2 = @$dbArRef2[0];
my $feat1 = $db1->gc->{$inputFeature};
my $feat2 = $db2->gc->{$inputFeature};
if ( !defined($feat1)  || !defined($feat2) ) 
{
  die "Feature $inputFeature does not exist in $dbfile1 or $dbfile2.";
}
my $cols1 = $db1->columns;
my $cols2 = $db2->columns;
if ( $cols1 != $cols2 ) 
{
  die "Number of cols in alignment1 = [$cols1] does not match alignment2 =[$cols2].";
}
print "$inputFeature-1=$feat1\n";
print "$inputFeature-2=$feat2\n";
my $matches = 0;
# TBD: make this a true input var
my @inputCat = ('H', 'E', 'L');
my @rocArray;
foreach my $cat (@inputCat)
{
  push @rocArray, Roc->new(cat=>$cat);
} 
for (my $col = 0; $col < $cols1; ++$col) {
  my $c1 = substr ($feat1, $col, 1);
  my $c2 = substr ($feat2, $col, 1);
  if ($c1 eq $c2) {
    ++$matches;
  }
  foreach my $rocRef (@rocArray)
  {
    $rocRef->evalColPred($c1, $c2);
  } 
}
foreach my $rocRef (@rocArray)
{
  $rocRef->toString();
} 
# Ignore "." columns in reference annotation
my $nullCols = $rocArray[0]->getNullCols();
my $effCols = $cols1 - $nullCols;
print "Orig cols=$cols1, Effective cols=$effCols\n";
print "Overall matches/cols = $matches/$effCols = ". $matches/$effCols*100 . "%\n";

sub help {
    print STDERR <<EOF;

$0: Compare column annotations between Stockholm alignments.
  The first annotation is assumed to be the reference.
  The second is assumed to be the prediction.
  Use "." to ignore those reference columns in calculations.

Usage: $0  <feature> <Stockholm alignment1> <Stockholm alignment2>
Options:
 
EOF
}
