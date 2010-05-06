#!/usr/bin/env perl -w

my $tmpseq;
END {
    if (defined($tmpseq) && -e $tmpseq) { unlink $tmpseq }
}

sub pathonly { my $path = shift; $path =~ s/(.*)\/([^\/]+)$/$1/; return $path }
sub nameonly { my $path = shift; $path =~ s/(.*)\/([^\/]+)$/$2/; return $path }
sub dosys { my $command = shift; warn $command, "\n"; system $command }

my $progname = nameonly ($0);
my $progpath = pathonly ($0);
my $dartpath = $ENV{'HOME'}."/dart";
my $dartbin = "$dartpath/bin";
my $dartperl = "$dartpath/perl";
my $dartsrc = "$dartpath/src";

my $handalign = "$dartperl/handalign.pl";
my $handcount = "$dartperl/handcount.pl";
my $tkfidem = "$dartbin/tkfidem";
my $hsmem = "$dartbin/hsmem";
my $inithsm = "$dartsrc/hsm/rndsubmat/HSM1a.max";

my $hsmopts = "-n";

my $usage = "\nUsage: $progname <seqfile> <output-dir>\n";

my (@log, $logfile);

my @argv;
while (@ARGV) {
    $arg = shift @ARGV;
    if ($arg =~ /^-/) {
	if ($arg eq "-log") { my $log; defined ($log = shift) or die $usage; push @log, $log }
	elsif ($arg eq "-logfile") { defined ($logfile = shift) or die $usage }
	else { die "$usage\nUnknown option: $arg\n\n" }
    } else { push @argv, $arg }
}
die $usage unless @argv == 2;

my ($seqfile, $outdir) = @argv;

my $tkfopts = join (" ", map ("-log $_", @log));
$tkfopts .= " -logfile $logfile" if defined $logfile;

mkdir $outdir unless -e $outdir;

my $bestloglike = -1e100;
my $insrate = .049505;
my $delrate = .050495;
for (my $iter = 1;; ++$iter) {
    my $align = "$outdir/align$iter.mul";
    my $tree = "$outdir/tree$iter.tree";
    my $binary = "$outdir/binary$iter.tree";
    my $distmat = "$outdir/distmat$iter.dist";
    my $index = "$outdir/index$iter.dist";
    my $idemfile = "$outdir/idem$iter";
    my $lasthsm = "$outdir/hsm".($iter-1);
    my $hsm = "$outdir/hsm$iter";
    if ($iter == 1) { dosys "cp $inithsm $lasthsm" }
    dosys "$handalign $tkfopts $seqfile -a $align -t $tree -m $distmat -ins $insrate -del $delrate -hsm $lasthsm";
    # get log-likelihood from alignfile
    my $loglike;
    local *ALIGN;
    open ALIGN, "<$align" or die "Opening $align: $!";
    while (<ALIGN>) {
	if (/Final\s*score\s*(\S+)\s*bits/) { $loglike = $1; last }
    }
    close ALIGN;
    die "No Final score in $align" if !defined $loglike;
    last if $loglike <= $bestloglike;
    $bestloglike = $loglike;
    # print index
    local *INDEX;
    open INDEX, ">$index" or die "Opening $index: $!";
    print INDEX "Alignment$iter\t$tree\t$align\n";
    close INDEX or die "Closing $index: $!";
    # get lambda, mu
    dosys "$tkfidem $tkfopts $index >$idemfile";
    local *IDEM;
    open IDEM, "<$idemfile" or die "Opening $idemfile: $!";
    while (<IDEM>) {
	if (/lambda:\s*(\S+)/) { $insrate = $1 }
	if (/mu:\s*(\S+)/) { $delrate = $1 }
    }
    close IDEM;
    # get submat
    dosys "$hsmem $hsmopts $tkfopts $lasthsm $index >$hsm";
}

print "lambda: $insrate\n";
print "mu:     $delrate\n";
