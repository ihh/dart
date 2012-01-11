#!/usr/bin/perl

use Getopt::Long;
use File::Basename;
use Carp;

use Stockholm;
use Stockholm::Database;

my ($progname) = fileparse($0);

my ($default_file_prefix, $fasta_suffix, $param_suffix, $tree_suffix) = qw(C1 P1.fastas p trees);
my $prefix = "./$default_file_prefix";

my $usage = "$progname -- convert a BaliPhy MCMC trace to a Stockholm alignment database\n";
$usage .= "\n";
$usage .= "Usage: $progname [-prefix <prefix>]\n";
$usage .= "\n";
$usage .= "Looks for '<prefix>.$fasta_suffix' (FASTA), '<prefix>.p' (params), '<prefix>.trees' (Newick)\n";
$usage .= "Default <prefix> is '$prefix'\n";
$usage .= "\n";

GetOptions ("prefix=s" => \$prefix) or die $usage;
$prefix .= "/$default_file_prefix" if -d $prefix;
my ($fasta_file, $param_file, $tree_file) = map ("$prefix.$_", $fasta_suffix, $param_suffix, $tree_suffix);

local *FASTA;
open FASTA, "<$fasta_file" or die "FASTA file '$fasta_file' not found";

warn "Parsing FASTA file '$fasta_file'\n";
local $_;
my $db = Stockholm::Database->new();
my (%seq, @seqorder, %iter2stock);
my ($iter, $name);
while (<FASTA>) {
    if (/^iterations = (\d+)/) {
	add ($db, $iter, \%seq, \@seqorder, \%iter2stock);
	$iter = $1;
	%seq = @seqorder = ();
	$name = undef;

	# parse the FASTA (it's sort of a shame this is duplicated from fasta2stockholm.pl ... FIXME?)
    } elsif (/^\s*>\s*(\S+)/) {
	$name = $1;
	die "Duplicate name: $name" if defined $seq{$name};
	push @seqorder, $name;

    } elsif (/\S/ && !defined $name) {
	warn "Ignoring: $_";

    } else {
	s/\s//g;
	$seq{$name} .= $_;
    }
}
add ($db, $iter, \%seq, \@seqorder, \%iter2stock);
close FASTA;

local *NEWICK;
open NEWICK, "<$tree_file" or die "Newick file '$tree_file' not found";

warn "Parsing Newick tree file '$tree_file'\n";
$iter = 0;
while (<NEWICK>) {
    if (exists $iter2stock{$iter}) {
	chomp;
	my $stock = $iter2stock{$iter};
	$stock->add_gf ('NH', $_);
    }
    ++$iter;
}
close NEWICK;

warn "Parsing parameter file '$param_file'\n";
local *PARAM;
open PARAM, "<$param_file" or die "Parameter file '$param_file' not found";
$_ = <PARAM>;
chomp;
my @param = split;
while (<PARAM>) {
    if (/\S/) {
	chomp;
	my @value = split;
	if (@value == @param) {
	    my %param2value = map (($param[$_], $value[$_]), 0..$#param);
	    my $iter = $param2value{'iter'};
	    if (exists $iter2stock{$iter}) {
		my $stock = $iter2stock{$iter};
		for my $i (0..$#param) {
		    $stock->add_gf ($param[$i], $value[$i]);
		}
	    }
	} else {
	    die "Failed to parse line of parameter file: $_";
	}
    }
}
close PARAM;

print $db->to_string;

sub add {
    my ($db, $iter, $seq_ref, $seqorder_ref, $iter2stock_ref) = @_;
    if (defined($iter) && @$seqorder_ref > 0) {
	my $stock = Stockholm->new();
	for my $seqname (@$seqorder_ref) {
	    $stock->add_row ($seqname, $seq_ref->{$seqname});
	}
	confess "Alignment not flush" unless $stock->is_flush;
	$stock->add_gf ('ID', $iter);
	$db->add_alignment ($stock);
	$iter2stock{$iter} = $stock;
    }
}
