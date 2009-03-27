#!/usr/bin/perl -w

my ($tmpseq, $tmptree, $tmpalign);
END {
    if (defined($tmpseq) && -e $tmpseq) { unlink $tmpseq }
    if (defined($tmptree) && -e $tmptree) { unlink $tmptree }
    if (defined($tmpalign) && -e $tmpalign) { unlink $tmpalign }
}

sub pathonly { my $path = shift; $path =~ s/(.*)\/([^\/]+)$/$1/; return $path }
sub nameonly { my $path = shift; $path =~ s/(.*)\/([^\/]+)$/$2/; return $path }
sub dosys { my $command = shift; warn $command, "\n"; system $command }

my $progname = nameonly ($0);
my $progpath = pathonly ($0);
my $dartpath = pathonly ($progpath);
my $dartbin = "$dartpath/bin";

my $tkfdistance = "$dartbin/tkfdistance";
my $weighbor = "$dartbin/weighbor";
my $tkfalign = "$dartbin/tkfalign";

my $usage = "\nUsage: $progname <seqfile>\n\t[-a <alignfile>] [-t <treefile>] [-m <distmatfile>]\n\t[-hsm <hsmfile>] [-ins <insrate>] [-del <delrate>]\n\t[-nowildcards] [-s <nsamples>] [-refine] [-keepinternal]\n\t[-log <loglevel>] [-logfile <logfile>]\n\n";

my ($distmat, $tree, $align, $hsmfile, $insrate, $delrate);
my (@log, $logfile);
my $alignopts = "";
my $keepinternal = 0;

my @argv;
while (@ARGV) {
    $arg = shift @ARGV;
    if ($arg =~ /^-/) {
	if ($arg eq "-m") { defined ($distmat = shift) or die $usage }
	elsif ($arg eq "-t") { defined ($tree = shift) or die $usage }
	elsif ($arg eq "-a") { defined ($align = shift) or die $usage }
	elsif ($arg eq "-hsm") { defined ($hsmfile = shift) or die $usage }
	elsif ($arg eq "-ins") { defined ($insrate = shift) or die $usage }
	elsif ($arg eq "-del") { defined ($delrate = shift) or die $usage }
	elsif ($arg eq "-log") { my $log; defined ($log = shift) or die $usage; push @log, $log }
	elsif ($arg eq "-logfile") { defined ($logfile = shift) or die $usage }
	elsif ($arg eq "-nowildcards") { $alignopts .= " -internal" }
	elsif ($arg eq "-s") { my $s; defined ($s = shift) or die $usage; $alignopts .= " -samples $s" }
	elsif ($arg eq "-refine") { $alignopts .= " -refine" }
	elsif ($arg eq "-keepinternal") { $keepinternal = 1 }
	else { die "$usage\nUnknown option: $arg\n\n" }
    } else { push @argv, $arg }
}
die $usage unless @argv == 1;

my ($seq) = @argv;
my $stub = $seq;
if (defined $align) { $stub = $align } else { $align = "$stub.mul" }
if (defined $tree) { $stub = $tree } else { $tree = "$stub.dnd" }
if (defined $distmat) { $stub = $distmat } else { $distmat = "$stub.dist" }
$tmpseq = "$seq.tmp";
$tmptree = "$tree.tmp";
$tmpalign = "$align.tmp";

unlink $distmat if -e $distmat;
unlink $tree if -e $tree;
unlink $align if -e $align;
unlink $tmpseq if -e $tmpseq;

local *SEQ;
local *TMP;
open SEQ, "<$seq" or die "Couldn't read $seq: $!";
open TMP, ">$tmpseq" or die "Couldn't write $tmpseq: $!";
while (my $tmpline = <SEQ>) {
    if ($tmpline =~ /^\s*>/) { $tmpline =~ tr/\/-/\@_/ }
    print TMP $tmpline;
}
close TMP;
close SEQ;

my $tkfopts = join (" ", map ("-log $_", @log));
$tkfopts .= " -logfile $logfile" if defined $logfile;
$tkfopts .= defined($hsmfile) ? " -hsm $hsmfile" : "";
$tkfopts .= " -birthrate $insrate" if defined $insrate;
$tkfopts .= " -deathrate $delrate" if defined $delrate;

# Get distance matrix
dosys "$tkfdistance $tkfopts $tmpseq >$distmat";

# Make tree by weighted neighbor-joining
dosys "$weighbor -i $distmat -o $tree";

# Check tree is properly terminated by a semicolon
my $newick = "";
local *TREE;
open TREE, "<$tree" or die $!;
while (<TREE>) { $newick .= $_ }
close TREE;
unless ($newick =~ /;\s*$/) {
    $newick =~ s/\s*$//;
    open TREE, ">$tree" or die $!;
    print TREE $newick, ";\n";
    close TREE or die $!;
}

# Run MCMC tree-alignment sampler
dosys "$tkfalign $tkfopts $alignopts $tree $tmpseq --save-tree $tmptree >$align";

# Update tree
rename $tmptree, $tree if -e $tmptree;

# discard internal sequences
unless ($keepinternal) {
    local *ALIGNIN;
    local *ALIGNOUT;
    open ALIGNIN, "<$align" or die "Couldn't read $align: $!";
    open ALIGNOUT, ">$tmpalign" or die "Couldn't write $tmpalign: $!";
    while (<ALIGNIN>) {
	print ALIGNOUT $_ unless /^root\s/ || /^\S+::\S+\s/;
    }
    close ALIGNOUT or die "Couldn't write $tmpalign: $!";
    close ALIGNIN;
    
    rename $tmpalign, $align if -e $tmpalign;
}
