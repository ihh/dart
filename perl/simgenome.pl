#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use strict;
use warnings;

use File::Spec;
use Carp;

use Stockholm;
use Newick;
use PhyloGram;

# get program name & directory
my $progname = $0;
$progname =~ s!.*/!!;

my $progdir = $0;
$progdir =~ s!/[^/]+$!!;

# get DARTDIR and grammar dir
my $dartdir = $ENV{'DARTDIR'};
if (!defined $dartdir) {
    if ($progdir =~ /\//) {
	$dartdir = $progdir;
	$dartdir =~ s!/[^/]+$!!;
    } else {
	$dartdir = "$progdir/..";
    }
}
my $grampath = "grammars";

# command-line args
my $treeFilename;
my $gramFilename = "simgenome.eg";
my $gsimroot = "trained-models/branch_model_singlet_2_12flymod100-pecan";
my $gsimbranch = "trained-models/branch_model_gotoh_mixture_2_12flymod30-pecan-95";
my $rndSeed;

# try looking in a few places for gsimulator
my $gsimdir = find_path ("gsimulator", ["/usr/local", "/usr/local/src", $ENV{'HOME'}], "-d", 0) || ".";

# usage text
my $usage;
$usage  = "Usage: $progname [options] [Newick tree file]\n";
$usage .= "Options:\n";
$usage .= "   -h         Print this help message\n";
$usage .= " -rnd         Seed random number generator to the current clock time\n";
$usage .= "-seed <int>   Seed random number generator to a particular value\n";
$usage .= "   -d <path>  Specify path to DART (default is \"$dartdir\")\n";
$usage .= "   -g <path>  Specify path to XRATE grammar (default is \"\$DARTDIR/$grampath/$gramFilename\")\n";
$usage .= "   -s <path>  Specify path to GSIMULATOR installation (default is \"$gsimdir\")\n";
$usage .= "  -rt <path>  Specify relative path to GSIMULATOR root transducer (default is \"$gsimroot\")\n";
$usage .= "  -bt <path>  Specify relative path to GSIMULATOR branch transducer (default is \"$gsimbranch\")\n";
$usage .= "\n";

# parse command-line opts
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-./) {
	if ($arg eq "-h") {
	    print $usage;
	    exit;
	} elsif ($arg eq "-rnd") {
	    $rndSeed = time;
	} elsif ($arg eq "-seed") {
	    defined ($rndSeed = shift) or die $usage;
	} elsif ($arg eq "-d") {
	    defined ($dartdir = shift) or die $usage;
	} elsif ($arg eq "-g") {
	    defined ($gramFilename = shift) or die $usage;
	} elsif ($arg eq "-s") {
	    defined ($gsimdir = shift) or die $usage;
	} elsif ($arg eq "-rt") {
	    defined ($gsimroot = shift) or die $usage;
	} elsif ($arg eq "-bt") {
	    defined ($gsimbranch = shift) or die $usage;
	} else {
	    die $usage, "Unknown option: $arg\n";
	}
    } else {
	die $usage if defined $treeFilename;
	$treeFilename = $arg;
    }
}

die $usage, "No tree specified\n" unless defined $treeFilename;

# print random number seed
my $rndSeedArg;
if (defined $rndSeed) {
    warn "[random number seed: $rndSeed]\n";
    $rndSeedArg = "-rndseed $rndSeed";
} else {
    warn
	"\n\n\n",
	"*** WARNING: no random number seed was supplied; default seed will be used. ***\n",
	"*** This means reruns of this program will produce EXACTLY THE SAME RESULT! ***\n",
	"\n\n\n\n";
    $rndSeedArg = "";
}

# find simgram
my $simgram = find_path ("simgram", "$dartdir/bin", "-x", 1);
my $simgram_args = "-gc #";   # set gap character to a little-used value, to suppress pointless warnings

# find grammar
my $gramfile_with_path = find_path ($gramFilename, "$dartdir/$grampath", "-e", 1);

# load tree
my $tree = Newick->from_file ($treeFilename);

# find gsimulator
my $gsimulator_with_path = find_path ("$gsimdir/bin/gsimulator.pl", [], "-x", 0);

# find gsimulator models
my $gsimroot_with_path = find_path ("$gsimdir/$gsimroot", [], "-e", 0);
my $gsimbranch_with_path = find_path ("$gsimdir/$gsimbranch", [], "-e", 0);

# load grammar
my $grammar = PhyloGram->from_file ($gramfile_with_path);

# find the tags that correspond to gsimulator
my (%gsim_args, %gsim_desc);
my $tags_sexpr = $grammar->grammar->meta->simgenome;
for my $sexpr ($tags_sexpr->find_all ("gsimulator")) {
    # check that we've got gsimulator & the models
    die "Can't find gsimulator directory '$gsimdir'\n" unless -d $gsimdir;
    my %gsim_files = ("gsimulator" => $gsimulator_with_path,
		      "gsimulator root model" => $gsimroot_with_path,
		      "gsimulator branch model" => $gsimbranch_with_path);
    while (my ($desc, $path) = each %gsim_files) {
	die "Can't find $desc -- are you sure gsimulator is in directory '$gsimdir' ?\n" unless defined($path) && length($path);
    }

    # add the gsimulator exec
    my $tag = $sexpr->row->value;
    my $args = $sexpr->find_all("args") ? join (" ", $sexpr->args->values) : "";
    $gsim_args{$tag} = "$gsimulator_with_path -r $gsimroot_with_path -b $gsimbranch_with_path $args";
    $gsim_desc{$tag} = "by lexicalized transducer simulation";
}

# find the tags that correspond to external simulation programs
for my $sexpr ($tags_sexpr->find_all ("external")) {
    my $tag = $sexpr->row->value;
    my @prog = $sexpr->executable->values;
    $gsim_args{$tag} = join (" ", @prog);
    $gsim_desc{$tag} = "using external program @prog";
}

# run simgram
my $simgram_cmd = "$simgram $simgram_args $rndSeedArg -g $gramfile_with_path -t $treeFilename";
warn "\nGenerating main template alignment and features by phylogrammar simulation (this may take a while...)\n $simgram_cmd\n";
local *SIMGRAM;
open SIMGRAM, "$simgram_cmd |";
my $stock = Stockholm->from_filehandle (\*SIMGRAM);
close SIMGRAM or die "While running simgram: $!\n";

# remove simgram's "#=GF observed-chain-counts" tags
delete $stock->gf->{"observed-chain-counts"};

# splice in gsimulator alignments
while (my ($gsim_tag, $gsim_exec) = each %gsim_args) {
    if (exists $stock->gc->{$gsim_tag}) {
	while (1) {
	    # look for flagged columns
	    my $annot = $stock->gc->{$gsim_tag};
	    last unless $annot =~ /([^\.]+)/g;

	    # find coords of flagged columns & flanking subalignments
	    my $match = $1;
	    my $rpos = pos($annot);
	    my $lpos = $rpos - length($match);

	    # get flanking subalignments
	    my $lflank = $stock->subalign (0, $lpos);
	    my $rflank = $stock->subalign ($rpos, $stock->columns - $rpos);

	    # run external simulator
	    my $gsim_cmd = "$gsim_exec $treeFilename";
	    warn "\nGenerating intergenic region at column ", $lpos+1, " $gsim_desc{$gsim_tag}\n $gsim_cmd\n";
	    local *GSIM;
	    open GSIM, "$gsim_cmd |";
	    my $gsim = Stockholm->from_filehandle (\*GSIM);
	    close GSIM or die "While running simgram: $!\n";

	    # remove all "#=GF" tags from externally simulated & right-flank alignments
	    $gsim->{'gf'} = {};
	    $rflank->{'gf'} = {};

	    # splice alignments together
	    $stock = $lflank;
	    $stock->concatenate ($gsim);
	    $stock->concatenate ($rflank);
	}

	# remove gsimulator tag
	delete $stock->gc->{$gsim_tag};
    }
}

# output alignment and exit
print $stock->to_string;
exit;

# subroutine to find a file, searching multiple paths
sub find_path {
    my ($file, $path, $file_test, $die_if_undef) = @_;
    $die_if_undef = 0 unless defined $die_if_undef;
#    carp "find_path($file,$path,$need_exec,$die_if_undef)";
    my @path = ref($path) eq 'ARRAY' ? @$path : ($path);
    if (defined $file) {
	for my $file_with_path ($file, map ("$_/$file", @path)) {
#	    warn "trying $file_with_path";
	    if (eval "$file_test \"$file_with_path\"") {
#		warn "find_path($file,$path,$need_exec,$die_if_undef): returning $file_with_path";
		return $file_with_path;
	    }
	}
    }
    die "Can't find $file in $path or elsewhere\n" if $die_if_undef;
#    carp "find_path($file,$path,$need_exec,$die_if_undef): returning undef";
    return undef;
}
