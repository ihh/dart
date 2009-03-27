#!/usr/bin/perl -w

my $tkfstalign = "tkfstalign";
my $grep = "grep";
my $fields = "fields";
my $GF_SC = "GF SC";   # Stockholm bitscore pattern
my $inf = 987654.321;

my ($tmin, $tmax, $tstep) = (0, 5, .1);
my ($lkmin, $lkmax, $lkstep) = (.1, .99, .1);
my ($lmmin, $lmmax, $lmstep) = (.01, 1, .01);
my ($skmin, $skmax, $skstep) = (.1, .99, .1);
my ($smmin, $smmax, $smstep) = (.01, 1, .01);
my ($spmin, $spmax, $spstep) = (.01, .99, .01);
my $verbose = 0;

my $usage = "\n";
$usage   .= "Usage: $0 <seqfile(s)>\n";
$usage   .= "   [-tmin  min_time]       [-tmax  max_time]       [-tstep  time_step]       (defaults $tmin, $tmax, $tstep)\n";
$usage   .= "   [-lkmin min_loop_kappa] [-lkmax max_loop_kappa] [-lkstep loop_kappa_step] (defaults $lkmin, $lkmax, $lkstep)\n";
$usage   .= "   [-lmmin min_loop_mu]    [-lmmax max_loop_mu]    [-lmstep loop_mu_step]    (defaults $lmmin, $lmmax, $lmstep)\n";
$usage   .= "   [-skmin min_stem_kappa] [-skmax max_stem_kappa] [-skstep stem_kappa_step] (defaults $skmin, $skmax, $skstep)\n";
$usage   .= "   [-smmin min_stem_mu]    [-smmax max_stem_mu]    [-smstep stem_mu_step]    (defaults $smmin, $smmax, $smstep)\n";
$usage   .= "   [-spmin min_stem_prob]  [-spmax max_stem_prob]  [-spstep stem_prob_step]  (defaults $spmin, $spmax, $spstep)\n";
$usage   .= "   [-t time] [-lk loop_kappa] [-lm loop_mu] [-sk stem_kappa] [-sm stem_mu] [-sp stem_prob]\n";
$usage   .= "   [-verbose]\n";
$usage   .= "\n";

my @argv;
while (@ARGV) {
    $arg = shift @ARGV;
    if ($arg =~ /^-/) {
	if ($arg eq "-tmin") { defined ($tmin = shift) or die $usage }
	elsif ($arg eq "-tmax") { defined ($tmax = shift) or die $usage }
	elsif ($arg eq "-tstep") { defined ($tstep = shift) or die $usage }
	elsif ($arg eq "-lkmin") { defined ($lkmin = shift) or die $usage }
	elsif ($arg eq "-lkmax") { defined ($lkmax = shift) or die $usage }
	elsif ($arg eq "-lkstep") { defined ($lkstep = shift) or die $usage }
	elsif ($arg eq "-lmmin") { defined ($lmmin = shift) or die $usage }
	elsif ($arg eq "-lmmax") { defined ($lmmax = shift) or die $usage }
	elsif ($arg eq "-lmstep") { defined ($lmstep = shift) or die $usage }
	elsif ($arg eq "-skmin") { defined ($skmin = shift) or die $usage }
	elsif ($arg eq "-skmax") { defined ($skmax = shift) or die $usage }
	elsif ($arg eq "-skstep") { defined ($skstep = shift) or die $usage }
	elsif ($arg eq "-smmin") { defined ($smmin = shift) or die $usage }
	elsif ($arg eq "-smmax") { defined ($smmax = shift) or die $usage }
	elsif ($arg eq "-smstep") { defined ($smstep = shift) or die $usage }
	elsif ($arg eq "-spmin") { defined ($spmin = shift) or die $usage }
	elsif ($arg eq "-spmax") { defined ($spmax = shift) or die $usage }
	elsif ($arg eq "-spstep") { defined ($spstep = shift) or die $usage }
	elsif ($arg eq "-t") { defined ($tmin = $tmax = shift) or die $usage }
	elsif ($arg eq "-lk") { defined ($lkmin = $lkmax = shift) or die $usage }
	elsif ($arg eq "-lm") { defined ($lmmin = $lmmax = shift) or die $usage }
	elsif ($arg eq "-sk") { defined ($skmin = $skmax = shift) or die $usage }
	elsif ($arg eq "-sm") { defined ($smmin = $smmax = shift) or die $usage }
	elsif ($arg eq "-sp") { defined ($spmin = $spmax = shift) or die $usage }
	elsif ($arg eq "-verbose") { $verbose = 1 }
	else { die "$usage\nUnknown option: $arg\n\n" }
    } else { push @argv, $arg }
}
die $usage unless @argv;

my @seqfile = @argv;
foreach my $seqfile (@seqfile) { die "Sequence file '$seqfile' not found" unless -e $seqfile }

my %tbest = map (($_ => undef), @seqfile);
my ($lk, $lm, $sk, $sm, $sp) = ($lkmax, $lmmax, $skmax, $smmax, $spmax);

iterate();


sub params {
    my $ll = $lk * $lm;
    my $sl = $sk * $sm;
    return "-ll $ll -lm $lm -sl $sl -sm $sm -sp $sp";
}

sub tscore {
    my ($seqfile, $t) = @_;
    my $params = params();
    my $sc = `$tkfstalign $params -t $t $seqfile | $grep "$GF_SC" | $fields 2`;
    chomp $sc;
    $sc = -$inf unless defined($sc) && $sc =~ /\d/;
    warn " [seqfile=$seqfile t=$t params=\"$params\" score=$sc]\n" if $verbose;
    return $sc;
}

sub tbest {
    my ($seqfile) = @_;
    my ($tbest, $scbest) = ($tmin, -$inf);
    for (my $t = $tmin; $t <= $tmax; $t += $tstep) {
	my $sc = tscore ($seqfile, $t);
	if ($sc > $scbest) { ($tbest, $scbest) = ($t, $sc) }
    }
    warn " Best t[$seqfile] = $tbest (score $scbest)\n";
    return $tbest;
}

sub optimise_times {
    foreach my $seqfile (@seqfile) {
	$tbest{$seqfile} = tbest ($seqfile);
    }
}

sub optimise_param {
    my ($name, $ref, $min, $max, $step) = @_;
    my ($best, $scbest) = ($min, -$inf);
    for (my $val = $min; $val <= $max; $val += $step) {
	$$ref = $val;
	my $sc = 0;
	foreach my $seqfile (@seqfile) {
	    $sc += tscore ($seqfile, $tbest{$seqfile});
	}
	if ($sc > $scbest) { ($best, $scbest) = ($val, $sc) }
    }
    $$ref = $best;
    warn " Best $name = $best (score $scbest)\n";
    return $scbest;
}

sub optimise_all {
    optimise_param ("lk", \$lk, $lkmin, $lkmax, $lkstep);
    optimise_param ("lm", \$lm, $lmmin, $lmmax, $lmstep);
    optimise_param ("sk", \$sk, $skmin, $skmax, $skstep);
    optimise_param ("sm", \$sm, $smmin, $smmax, $smstep);
    my $sc = optimise_param ("sp", \$sp, $spmin, $spmax, $spstep);
    return $sc;
}

sub iterate {
    my $last_sc = -$inf;
    my $round = 0;
    while (1) {
	warn "Starting round ", ++$round, "\n";
	optimise_times();
	my $sc = optimise_all();
	last if $sc <= $last_sc;
	$last_sc = $sc;
    }
    print params(), "\n";
}
