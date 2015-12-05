#!/usr/bin/env perl

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use strict;
use warnings;

use Getopt::Long qw(GetOptions);

use Stockholm;

my ($infile, $outfile, $tree, $dna);

my $dartdir = $ENV{'DARTDIR'} || ($FindBin::Bin . "/..");
my $muscle = `which muscle` || "$dartdir/bin/muscle";
chomp $muscle;

my $usage = "\n"
    ."Usage: protpaleo.pl <infile> [more options...]\n"
    ."\n"
    ."   <infile>             (FASTA or Stockholm format)\n"
    ."   -dna                 (if input is DNA)\n"
    ."   -tree <tree string>  (New Hampshire format)\n"
    ."   -out <outfile>\n"
    ."\n";

@ARGV || die $usage;
unshift @ARGV, "-in";
GetOptions ("in=s" => \$infile,
	    "out=s" => \$outfile,
	    "tree=s" => \$tree,
	    "dartdir=s" => \$dartdir,
	    "muscle=s" => \$muscle,
	    "dna" => \$dna)
    or die $usage;

defined($infile) or die $usage;

my $dartbin = "$dartdir/bin";

my $xrate = "$dartbin/xrate";
my $protpal = "$dartbin/protpal";

my $nullgram = "$dartdir/grammars/" . ($dna ? "nullrna.eg" : "nullprot.eg");

my $muscle_opts = $dna ? "-maxiters 1 -diags" : "-maxiters 1 -diags -sv -distance1 kbit20_3";

# Fasta or Stockholm?
my $got_stock = 0;
local *IN;
open IN, "<$infile" or die "Couldn't open input file $infile: $!";
while (<IN>) {
    if (/^\s*# STOCKHOLM"/) {
	$got_stock = 1;
	last;
    } elsif (/\S/ && /^\s*[^#]/) {
	last;
    }
}
close IN;

# Get guide alignment
my $alignfile;
if ($got_stock) {
    $alignfile = $infile;
} else {
    $alignfile = "$infile.align.stockholm";
    my $alignfile_fasta = "$infile.align.fasta";
    sysrun ("$muscle $muscle_opts -in $infile -out $alignfile_fasta");
    fasta2stockholm ($alignfile_fasta, $alignfile);
}

# Get tree
$tree = $tree || stockholm_tree($alignfile);
if (!defined $tree) {
    my $treealignfile = "$alignfile.treealign.stockholm";
    sysrun ("$xrate -g $nullgram -noa $alignfile >$treealignfile");
    $tree = stockholm_tree ($treealignfile);
}

# Write tree
my $treefile = "$infile.newick";
local *TREE;
open TREE, ">$treefile";
print TREE $tree;
close TREE;

# Run protpal
my $infile_opts = $got_stock ? "-stk $infile" : "-fa $infile";
my $outfile_opts = defined($outfile) ? ">$outfile" : "";
my $dna_opts = $dna ? "-dna" : "";
sysrun ("$protpal -ga $alignfile -tf $treefile -g $nullgram $dna_opts -eri $outfile_opts");


# Helper functions
sub fasta2stockholm {
    my ($fasta, $stockholm) = @_;
    my %seq;
    my @name;
    my $name;
    local *FASTA;
    open FASTA, "<$fasta" or die "Couldn't open FASTA file '$fasta': $!";
    while (<FASTA>) {
	if (/^\s*>\s*(\S+)/) {
	    $name = $1;
	    die "Duplicate name: $name" if defined $seq{$name};
	    push @name, $name;
	} else {
	    if (/\S/ && !defined $name) {
		warn "Ignoring: $_";
	    } else {
		s/\s//g;
		$seq{$name} .= $_;
	    }
	}
    }
    close FASTA;

    my $length;
    my $lname;
    foreach my $name (@name) {
	my $l = length $seq{$name};
	if (defined $length) {
	    die "Sequences not all same length ($lname is $length, $name is $l)" unless $length == $l;
	} else {
	    $length = length $seq{$name};
	    $lname = $name;
	}
    }

    local *STOCK;
    open STOCK, ">$stockholm" or die "Couldn't write to file $stockholm: $!";
    print STOCK "# STOCKHOLM 1.0\n";
    foreach my $name (@name) {
	print STOCK $name, " ", $seq{$name}, "\n";
    }
    print STOCK "//\n";
    close STOCK;
}

sub stockholm_tree {
    my ($stockfile) = @_;
    my $stock = Stockholm->from_file ($stockfile);
    return @{$stock->gf_NH} ? $stock->NH : undef;
}

sub sysrun {
    my ($cmd) = @_;
    warn "Running $cmd\n";
    system $cmd;
}
