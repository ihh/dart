#!/usr/local/bin/perl -w

my $tmpseq;
END {
    if (defined($tmpseq) && -e $tmpseq) { unlink $tmpseq }
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
my $tkfidem = "$dartbin/tkfidem";

my $usage = "\nUsage: $progname <alignment#1> [<alignment#2> ...]\n\t[-log <loglevel>] [-logfile <logfile>]\n\n";

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
die $usage unless @argv > 0;

my @align = @argv;

my $tkfopts = join (" ", map ("-log $_", @log));
$tkfopts .= " -logfile $logfile" if defined $logfile;

my $index = join ("-", @align);
$index =~ s/\.mul//g;
$index .= ".idx";
open INDEX, ">$index";
foreach my $align (@align) {
    my $dist = "$align.dist";
    my $tree = "$align.tree";

    dosys "$tkfdistance $tkfopts $align >$dist";
    dosys "$weighbor -i $dist -o $tree";

    print INDEX "$align\t$tree\t$align\n";
}
close INDEX;

dosys "$tkfidem $tkfopts $index";
