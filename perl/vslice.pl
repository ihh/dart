#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm;

# command-line params
my $refseqname;

# usage
my $usage = "Usage: $0 <Stockholm alignment>
            [-h|--help]          Print this message.
            [-r|--refseq <NAME>] Treat <NAME> as the reference sequence, instead of the first one.
\n";

my @argv;
while (@ARGV) {
  my $arg = shift;
  if ($arg eq "-h" || $arg eq "--help") { print $usage; exit; }
  elsif ($arg eq "-r" || $arg eq "--refseq") { defined ($refseqname = shift) or die $usage }
  else { push @argv, $arg }
}

unless (@argv) {
  @argv = ('-');
  warn "[waiting for alignment on standard input]\n";
}

my ($file) = shift @argv;
warn "ignoring files: @argv\n" if @argv;

my $stock = Stockholm->from_file ($file);

if ($stock->sequences) {
    $refseqname = $stock->seqname->[0] unless defined $refseqname;
    my $refseq = $stock->seqdata->{$refseqname};
    die "Can't find reference sequence '$refseqname'\n" unless defined $refseq;
    my $gapCharRe = Stockholm::gapCharsRegexp();
    my @gappy_cols = grep (substr ($refseq, $_, 1) =~ /[$gapCharRe]/, 0..$stock->columns-1);
    $stock->drop_columns (@gappy_cols);
}

print $stock->to_string;

