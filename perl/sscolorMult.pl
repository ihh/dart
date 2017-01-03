#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

=head1 NAME

sscolorMult.pl - call sscolor.pl repeatedly to generate multiple interlinked HTML files

=head1 SYNOPSIS

 sscolorMult.pl -h

 sscolorMult.pl [-o PREFIX] FILENAME

=head1 DESCRIPTION

See documentation for sscolor.pl.

=head1 OPTIONS

=over 12

=item C<-h>

Print a short help message.

=item C<-o>

Select prefix for output files (passed to sscolor.pl as the C<-out> option).

=back

=cut

use strict;
use Getopt::Long;
use Stockholm;

# Read command line arguments

my $out = "out";
my ($help, @paths);

GetOptions( "h"     => \$help,
            "o=s"   => \$out
          );
	     
my $file = shift;

if( $help or not $file ) {
    help();
    exit(1);
}

# Main 

my $stk = Stockholm->from_file($file);
my $numSeq = $stk->sequences();
for (my $i=1; $i<= $numSeq; $i++)
{
  system("sscolor.pl -link -ref $i -out $out $file > $out.$i.html");
}

# Subroutines

sub help {
    print STDERR <<EOF;

$0: 
Calls sscolor.pl multiple times to output a linked set of html colorings
using every sequence in the MSA as the reference sequence.

Usage: $0 [options] stockholmAlignment
    Options
        -h             : show this help
	-o <fileName>  : output base filename [default='$out']
EOF
}
