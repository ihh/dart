#!/usr/bin/perl -w

my $stockSuffix = ".stockholm";
my $xrangeSuffix = ".xrange";
my $xrateSuffix = ".xrate";

my $bindir = "../../../bin";
my $xrate = "$bindir/irrev_xrate";

local $DIR;
opendir DIR, '.';
my @dir = grep -f $_, readdir DIR;
closedir DIR;

my %dir = map (($_ => 1), @dir);
my @prefix = map (/^(\S+)$stockSuffix$/ && $dir{"$1$xrangeSuffix"} ? ($1) : (), @dir);

my $ok = 1;
foreach my $prefix (@prefix) {
    my ($stockFile, $xrangeFile, $xrateFile) = map ("$prefix$_", $stockSuffix, $xrangeSuffix, $xrateSuffix);

    local *STOCK;
    open STOCK, "<$stockFile" or die "'$stockFile': $!";
    my $opts = "";
    my $fail = 0;
    while (my $stock = <STOCK>) {
	if ($stock =~ /^\s*\#=GF\s+CC\s+xrate\s+(.*)$/) {
	    $opts = $1;
	} elsif ($stock =~ /^\s*\#=GF\s+CC\s+fail$/) {
	    $fail = 1;
	}
    }
    close STOCK;

    my $command = "$xrate $opts $stockFile >$xrateFile";
    warn "[Comparing '$command' with '$xrangeFile']\n";
    system $command;

    local *XRATE;
    open XRATE, "<$xrateFile" or die "'$xrateFile': $!";
    my @xrateOutput = <XRATE>;
    close XRATE;

    local *XRANGE;
    open XRANGE, "<$xrangeFile" or die "'$xrangeFile': $!";
    my @desiredOutput = <XRANGE>;
    close XRANGE;

    my $testok = 1;
    if (@desiredOutput != @xrateOutput) {
	$testok = 0;
    } else {
	for (my $i = 0; $i < @desiredOutput; ++$i) {
	    my $x = $xrateOutput[$i];
	    my $d = $desiredOutput[$i];
	    $x =~ s/^\s+(.*?)\s+$/$1/;
	    $d =~ s/^\s+(.*?)\s+$/$1/;
	    my @x = split /\s+/, $x;
	    my @d = split /\s+/, $d;
	    if (@x != @d) {
		warn "   [Line $i of xrate output has ", @x+0, " fields; expected ", @d+0, "]\n";
		$testok = 0;
	    } else {
		for (my $j = 0; $j < @d; ++$j) {
		    my $dj = $d[$j];
		    if ($dj =~ /^([0-9\.\-]+)\:([0-9\.\-]+)$/) {
			my ($min, $max) = ($1, $2);
			if ($min > $max) { ($min, $max) = ($max, $min) }
			unless ($x[$j] >= $min && $x[$j] <= $max) {
			    warn "   [Field $j of line $i is $x[$j]; permissible range is $min to $max]\n";
			    $testok = 0;
			}
		    } else {
			if ($x[$j] ne $dj) {
			    warn "   [Field $j of line $i is $x[$j]; should be $dj]\n";
			    $testok = 0;
			}
		    }
		}
	    }
	}
    }

    if ($testok && !$fail) {
	warn "  [comparison ok]\n";
    } elsif ($testok && $fail) {
	warn "  [comparison ok, EVEN THOUGH FAILURE WAS EXPECTED -- update file '$stockFile']\n";
    } elsif (!$testok && $fail) {
	warn "  [test failed, but failure was expected]\n";
    } else {  # !$testok && !$fail
	warn "Desired xrate output from file '$xrangeFile'\n", "<BEGIN>\n", @desiredOutput, "<END>\n";
	warn "Actual xrate output from command '$command'\n", "<BEGIN>\n", @xrateOutput, "<END>\n";
	$ok = 0;
    }
}

print $ok ? "ok\n" : "not ok\n";
