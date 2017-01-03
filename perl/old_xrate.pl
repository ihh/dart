#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

my $hsm2xgram = "hsm2xgram.pl";
my $xgram2hsm = "xgram2hsm.pl";
my $xrate = "xrate";

my $model = "-p rev";

my $tmpIn = "tmp-input.eg";
my $tmpOut = "tmp-output.eg";

my $usage = "$0: wrapper for XRATE preserving old HSM file format\n";

my @argv;
while (my $arg = shift) {
    if ($arg eq "-i" || $arg eq "-init" || $arg eq "--initialise") {
	my $hsmFile = shift;
	system "$hsm2xgram $hsmFile >$tmpIn";
	$model = "-g $tmpIn";
    } elsif ($arg eq "-irrev") {
	$model = "-p irrev";
    } elsif ($arg eq "-c") {
	my $classes = shift;
	$classes = "" if $classes == 1;
	$model = "-p aa$classes";
    } else {
	push @argv, $arg;
    }
}

system "$xrate $model --noannotate --train $tmpOut @argv >/dev/null";
system "$xgram2hsm $tmpOut";

