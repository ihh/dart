#!/usr/bin/env perl -w

use strict;
use CGI qw/:standard *table *pre *span/;

my $PS = "PS";
my $SS = "SS";
my $SC = "SC";
my $PS_CONS = "PS_cons";
my $SS_CONS = "SS_cons";
my $SN = "SN";

my $ps_cons_hdr = "#=GC $PS_CONS";
my $ss_cons_hdr = "#=GC $SS_CONS";
my $sn_hdr = "#=GC $SN";
my $sc_hdr = "#=GF $SC";

# Basepair compatibility levels
my $BP_COMPAT_NO_SUB = 0;
my $BP_COMPAT_SGL_SUB = 1;
my $BP_COMPAT_DBL_SUB = 2;
my $BP_NON_COMPAT = 3;

# Create hash maps for secondary structure characters.
# Left to right character
my %rchar = ('('=>')','<'=>'>','['=>']');
# Right to left character
my %lchar = map (($rchar{$_}=>$_), keys(%rchar));
# Gap characters
my %gapchar = map (($_=>1), '-', '.', '_', ',');

my $refSeqNum = 1;
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/; # Remove path
my $verbose = 0;
my $makeLinks = 0;
my $outFile = "out";

my $usage = "$progname:  Color the secondary structure of alignments in Stockholm file.\n";
$usage .= "   Secondary structure should be in a line prefixed by '#=GC SS_cons'.\n";
$usage .= "   Substitutions are regarded with respect to the reference sequence.\n";
$usage .= "   Coloring of nucleotides is as follows:\n";
$usage .= "       Green = compatible Watson-Crick basepair with double substitution\n";
$usage .= "       Blue = compatible W-C basepair with a single substitution\n";
$usage .= "       Red = substitution resulting in a non-compatible basepair\n";
$usage .= "       Black = no substitution\n";
$usage .= "       Grey = non-basepair\n";
$usage .= "       Light Purple = non-basepair with a substitution\n";
$usage .= "\n";
$usage .= "Usage: $progname <Stockholm file>\n";
$usage .=   "             [-h] print this message\n";
$usage .=   "      [-ref <n> ] reference sequence number (default = 1)\n";
$usage .=   "          [-link] add links so viewer can select reference sequence (use with sscolorMult.pl)\n";
$usage .=   "   [-out <file> ] output file basename used in links (default = 'out')\n";
$usage .=   "\n";
$usage .=   "For more help, type the following:\n";
$usage .=   "perldoc $0\n";
$usage .=   "\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") 
    {
      die $usage;
    } 
    elsif ($arg eq "-v") 
    {
      $verbose = 1;
    } 
    elsif ($arg eq "-ref") 
    {
      defined ($refSeqNum = shift) or die $usage;
    } 
    elsif ($arg eq "-link") 
    {
      $makeLinks = 1;
    } 
    elsif ($arg eq "-out") 
    {
      defined ($outFile = shift) or die $usage;
    } 
    else 
    {
      push @argv, $arg;
    }
}
@argv = ('-') unless @argv;
die $usage unless @argv == 1;
my ($stockfile) = @argv;

# Define the document style sheet and start printing out the html.
my $style=<<END;
    .noPair {
      color: #999999; /* med gray */
    }
    .noPairSub {
      color: #9900CC; /* dark purple */
    }
    .compatSgl {
      color: #0000FF; /* blue */
    }
    .compatDbl {
      color: #00FF00; /* green */
    }
    .nonCompat {
      color: #FF0000; /* red */
    }
    .noSub {
      color: #000000; /* black */
    }
END
print start_html(-title=>$stockfile, -style=>{-code=>$style}),"\n";

local *STOCK;
open STOCK, "<$stockfile" or die "Couldn't open '$stockfile': $!";
while (1) 
{
  my (@seqname, %seq, %ss, %gc, %gr, $score);
  my ($ps_cons, $ss_cons, $preamble);
  my $seqlen = 0;
  my $found_separator = 0;

  # Read in next alignment in Stockholm file.
  while (<STOCK>) 
  {
    # "#=GR...SS" by-seq secondary structure line
    if (/^\s*\#=GR\s*(\S+)\s*$SS\s*(\S+)\s*$/) 
    { 
      $ss{$1} = "" unless exists $ss{$1}; 
      $ss{$1} .= $2; 
    }
    # miscellaneous "#=GR" lines
    elsif (/^\s*\#=GR\s*(\S+)\s*(\S+)\s*(\S+)\s*$/) 
    { 
      $gr{$2} = {} unless exists $gr{$1}; 
      # %gr will store a hash reference.
      $gr{$1}->{$2} = "" unless exists $gr{$1}->{$2}; 
      $gr{$1}->{$2} .= $3 
    }  
    # "#=GC...PS_cons" consensus primary sequence
    elsif (/^\s*\#=GC\s*$PS_CONS\s*(\S+)\s*$/) 
    { 
      $ps_cons = "" unless defined $ps_cons; 
      $ps_cons .= $1; 
    }  
    # "#=GC...SS_cons" consensus secondary structure
    elsif (/^\s*\#=GC\s*$SS_CONS\s*(\S+)\s*$/) 
    { 
      $ss_cons = "" unless defined $ss_cons; 
      $ss_cons .= $1 
    }  
    # miscellaneous "#=GC" lines
    elsif (/^\s*\#=GC\s*(\S+)\s*(\S+)\s*$/) 
    { 
      $gc{$1} = "" unless exists $gc{$1}; 
      $gc{$1} .= $2; 
    }  
    # score
    elsif (/^\s*\#=GF\s*$SC\s*(\S+)\s*$/) 
    { 
      $score = $1; 
    }  
    # unrecognised line starting with '#'; append to preamble
    elsif (/^\s*\#/) 
    { 
      $preamble = "" unless defined $preamble; 
      $preamble .= $_ 
    }  
    # Found alignment separator; exit loop
    elsif (/^\s*\/\//) 
    { 
      $found_separator = 1; 
      last; 
    }  
    elsif (/^\s*(\S+)\s*(\S+)\s*$/) 
    { 
      unless (exists $seq{$1}) 
      { 
	$seq{$1} = ""; 
	push @seqname, $1 
      } 
      $seq{$1} .= $2; 
      $seqlen = length($seq{$1}) 
    }
    elsif (/\S/) 
    { 
      warn "Ignoring line: $_" 
    }
  } # while (<STOCK>)...

  if (@seqname) 
  {
    if (defined($ss_cons) && length($ss_cons)) 
    {
      my @pairs;
      my $columns = length($ss_cons);
      my $rows = @seqname;
      my @lpos;
      # Get coords of all basepairs.
      for (my $i = 0; $i < length($ss_cons); ++$i) 
      {
	my $c = substr ($ss_cons, $i, 1);
	# Found left character.
	if (exists $rchar{$c}) 
	{  
	  push @lpos, $i;
	} 
	# Found right character
	elsif (exists $lchar{$c}) 
	{  
	  unless (@lpos) 
	  { 
	    die "Too many >'s in $SS_CONS string: $ss_cons\n" 
	  }
	  my $l = pop @lpos;
	  my $r = $i;
#	  warn "l=$l r=$r\n";
	  push @pairs, {coord => [$l,$r], level => $BP_COMPAT_NO_SUB};
	}
      }
      if (@lpos) 
      { 
	die "Too many <'s in $SS_CONS string: $ss_cons\n" 
      }
      # Determine maximum sequence name length
      my $maxNameLen = max(map(length($_), $SS_CONS, @seqname));
      my $fmt = join("", "%-", $maxNameLen, "s ");
      # Set reference seq array
      if ($refSeqNum > $rows ) 
      { 
	die "Reference sequence number [$refSeqNum] is greater than number of sequences [$rows]\n";
      }
      my @refSeq = split("", $seq{$seqname[$refSeqNum-1]}); 
      # Print alignment
      print start_pre(),"\n";
      for (my $row = 0; $row < $rows; ++$row) 
      {
	my $isRefRow = ($row == 0 ? 1 : 0);
	my $seqName = $seqname[$row];
	print $makeLinks ? a({-href=>$outFile . "." . ($row+1) . ".html"}, $seqName) : $seqName;
	printf(" " x ($maxNameLen-length($seqName) + 1) );
	my $rowseq = $seq{$seqName};
	# Loop thru each character in a row.
	for (my $i=0; $i < length($rowseq); ++$i)
	{
	  my $qChar = substr($rowseq, $i, 1);
	  my $refChar = $refSeq[$i];
	  my $class = "noPair";
	  my $level;
	  # Determine whether the character is a base pair member and
	  #   the level of base pair compatibility.
	  # $posL, $posR refer to the left, right positions of a basepair.
	  my ($posLflag, $posRflag) = (0,0);
	  foreach my $pair (@pairs)
	  {
	    my $posL = $pair->{coord}[0];
	    my $posR = $pair->{coord}[1];
	    if ($i == $posL)
	    {
	      $posLflag = 1;
	      if ($isRefRow)
	      {
		$level = $BP_COMPAT_NO_SUB;
	      }
	      else
	      {
		print "[L:$qChar]" if $verbose;
		my $qChar2 = substr($rowseq, $posR,1);
		my $refChar2 = $refSeq[$posR];
		$level = getBpCompat($qChar, $qChar2, $refChar, $refChar2);
		print ("[$qChar$qChar2$refChar$refChar2:$level]") if $verbose;
	      }
	      $pair->{level} = $level;
	      last;
	    }
	    elsif ($i == $posR)
	    {
	      $posRflag = 1;
	      $level = $pair->{level};
	      print "[i=$i, posL=$posL, posR=$posR]" if $verbose;
	      print "[R:$qChar,$level]" if $verbose;
	      last;
	    }
	  }
	  # Character is a pair member.
	  if ( $posLflag || $posRflag)
	  {
	    if ($level == $BP_COMPAT_NO_SUB)
	    {
	      $class='noSub';
	    }
	    elsif ($level == $BP_COMPAT_SGL_SUB)
	    {
	      $class='compatSgl';
	    }
	    elsif ($level == $BP_COMPAT_DBL_SUB)
	    {
	      $class='compatDbl';
	    }
	    elsif ($level == $BP_NON_COMPAT)
	    {
	      $class='nonCompat';
	    }
	  }
	  else
	  {
	    if ($qChar ne $refChar)
	    {
	      $class ='noPairSub';
	    }
	  }
	  print span({-class=>$class}, $qChar);
	}
	print "\n";
      }
      printf(join("",$fmt,"%s"), $SS_CONS, $ss_cons);
      print "\n", end_pre(),"\n";
    } # end of block that looks at SS_cons line
  } # if (@seqname) ...

  # if no alignment separator, then quit the loop
  last unless $found_separator;
}
close STOCK;
print end_html();

# Determine the compatibility of a query sequence base pair relative to 
#   a reference sequence base pair.
# Output: the level of compatibility.
sub getBpCompat
{
  # q1, q2 are the query seq bases.
  # r1, r2 are the reference seq bases.
  my ($q1, $q2, $r1, $r2) = @_; 
  my $level = $BP_COMPAT_NO_SUB;
  my $numSub = 0;
  if (uc($q1) ne uc($r1))
  {
    $numSub++;
  }
  if (uc($q2) ne uc($r2))
  {
    $numSub++;
  }
  if ($numSub > 0)
  {
    if (isValidPair($q1, $q2))
    {
      if ($numSub == 1)
      {
	$level = $BP_COMPAT_SGL_SUB;
      }
      else
      {
	$level = $BP_COMPAT_DBL_SUB;
      }
    }
    else
    {
      $level = $BP_NON_COMPAT;
    }
  }

  return $level;
}

# Determines whether a DNA/RNA basepair is valid.
sub isValidPair
{
  my ($b1, $b2) = @_;
  my %validPairs = ( A_T => 1, A_U => 1, G_C => 1, G_U => 1); 
  my $isValid = 0;
  my $pair1 = $b1."_".$b2;
  my $pair2 = $b2."_".$b1;
  if (exists($validPairs{uc($pair1)}) || exists($validPairs{uc($pair2)}))
  {
    $isValid = 1;
  }
  return $isValid;
}

# Determine the list maximum value.
sub max 
{
  my ($x, @y) = @_;
  foreach my $y (@y) 
  { 
    $x = $y if $y > $x 
  }
  return $x;
}


=head1 NAME

sscolor.pl - colorize RNA multiple alignments using HTML

=head1 SYNOPSIS

 sscolor.pl -h

 sscolor.pl FILENAME [-ref N] [-out PREFIX]

 cat FILENAME | sscolor.pl [options]

=head1 DESCRIPTION

The program uses HTML colors to highlight compensatory mutations in a multiple alignment of noncoding RNA.

The input alignment must be in Stockholm format.
The output is in HTML.

For the purposes of this program,
a position is defined to display a compensatory mutation if (i) it is base-paired with another position;
(ii) the pair of nucleotides at those two positions comprises a Watson-Crick or wobble base-pair;
(iii) at least one of the two nucleotides is different from the homologous nucleotide in a reference sequence.

=head1 INPUT FORMAT

The input file(s) must be in Stockholm format, defined here:

L<http://biowiki.org/StockholmFormat>

A list of utilities for format conversions and other operations on Stockholm files can be found here:

L<http://biowiki.org/StockholmTools>

=head1 OUTPUT FORMAT

Output in HTML is printed on the standard output.

The C<-link> and C<-out> options can be used to create a family of HTML pages,
with each one using a different reference sequence.
Use the accompanying sscolorMult.pl script to do this automatically.

Coloring of nucleotides is as follows:

=over 2

=item Green

Consensus structure implies Watson-Crick or wobble base-pair.

Double substitution relative to reference sequence.

=item Blue

Consensus structure implies Watson-Crick or wobble base-pair.

Single substitution relative to reference sequence.

=item Black

Consensus structure implies Watson-Crick or wobble base-pair.

No substitution relative to reference sequence.

=item Red

Consensus structure implies non-canonical base-pair.

=item Grey

Consensus structure implies position is single-stranded.

No substitution relative to reference sequence.

=item Light purple

Consensus structure implies position is single-stranded.

Single substitution relative to reference sequence.

=back

=head1 OPTIONS

=over 12

=item C<-h>

Prints a short help message.

=item C<-ref N>

Allows the user to select which row number will be used for the reference sequence.

=item C<-link>

Creates internal links for each sequence name, which the user can click on (in the web browser) to select different reference sequences.

This option presupposes that a separate HTML file has been created for each row of the alignment, wherein that row has been selected as the reference sequence.
To generate these files automatically, use the sscolorMult.pl script.

=item C<-out PREFIX>

Sets the filename prefix for the internal links generated by the C<-link> option.
Used by the sscolorMult.pl script.

=back

=head1 AUTHOR

Yuri Bendana - <ybendana@berkeley.edu>

=head1 LICENSE

This code is released under the GNU Public License v3.0.

=cut
