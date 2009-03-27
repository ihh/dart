#!/usr/bin/perl -w

# script to simulate pairwise alignments using handalign's model, and measure max deviation from diagonal

my $prog = $0;
$prog =~ s/.*\///;

my $reps = 100;  # number of trials
my $dmin = .1;  # minimum deletion rate
my $dmax = 1;  # maximum deletion rate
my $dstep = .1;  # step size for deletion rate
my $gmin = 1;  # minimum mean gap length
my $gmax = 101;  # maximum mean gap length
my $gstep = 10;  # step size for mean gap length
my $kappa = .99;  # insertion/deletion rate ratio; equal to 1-1/(L+1) where L = mean length
my $anclen;  # ancestral sequence length (if defined then ancestral seq.len. is fixed at this value, rather than geometrically distributed)

my $showpath = 0;
my $cts = 0;
my $idcount = 0;

my @params = qw(reps dmin dmax dstep gmin gmax gstep kappa);
my $regexp = "^-(" . join ("|", @params) . ")\$";
my $defaults = join (", ", map ("$_(".eval("\$$_").")", @params));

my @outfields = qw(deletion_rate mean_gap_size delete_open_prob insert_open_prob delete_extend_prob insert_extend_prob longest_seq_length desc_len_minus_anc_len max_dev max_skew);

my $usage = "\n";
$usage .= "Usage: $prog [-param_name param_val]\n";
$usage .= "\n";
$usage .= "Parameters (& defaults):\n";
$usage .= " $defaults\n";
$usage .= "\n";
$usage .= "Shortcuts & other options:\n";
$usage .= " -cts         (do exact continuous-time simulation instead of a transducer approximation)\n";
$usage .= "  -id          (track insertion & deletion event counts when in continuous-time mode; useful for checking detailed balance)\n";
$usage .= " -path        (show run-length encoded state path as an extra " . (@outfields+1) . "th field)\n";
$usage .= " -d d_value   (sets dmin=dmax=d_value)\n";
$usage .= " -g g_value   (sets gmin=gmax=g_value. NB this is the mean deletion size, not the mean insertion size)\n";
$usage .= " -l length    (fixes ancestral sequence length)\n";
$usage .= "\n";
$usage .= "Output consists of " . (@outfields+0) . " tab-delimited fields:\n";
$usage .= " " . join ("  ", map ("$outfields[$_]($_)", 0..$#outfields)) . "\n";
$usage .= "\n";
$usage .= "Here max_dev = max |y-diag(x)|, where diag(x) = x * L_y / L_x;  whereas max_skew = max |y-x|.\n";
$usage .= "Note that some of the fields are related:\n";
$usage .= " delete_extend_prob = 1 - 1 / (1 + mean_gap_size)        [so mean_gap_size param is actually *one less than* the true mean gap size]\n";
$usage .= " insert_extend_prob = kappa * delete_extend_prob\n";
$usage .= "   delete_open_prob = 1 - exp(-deletion_rate)\n";
$usage .= "   insert_open_prob = 1 - exp(-insertion_rate)\n";
$usage .= "     insertion_rate = deletion_rate * (1 - delete_extend_prob) / (1 - insert_extend_prob)\n";
$usage .= "\n";
$usage .= "Example usage:\n";
$usage .= " ./simband.pl -d .1 -g 1 -reps 10000 -l 200 | fields 7 8 | xybins.pl 2 | xgraph\n";
$usage .= "\n";

# parse command-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-./) {
	if ($arg eq "-h") {
	    print $usage;
	    exit;
	} elsif ($arg eq "-d") {
	    defined ($dmin = shift) or die $usage;
	    $dmax = $dmin;
	} elsif ($arg eq "-g") {
	    defined ($gmin = shift) or die $usage;
	    $gmax = $gmin;
	} elsif ($arg eq "-l") {
	    defined ($anclen = shift) or die $usage;
	} elsif ($arg =~ /$regexp/) {
	    my $val;
	    defined ($val = shift) or die $usage;
	    eval ("\$$1 = $val");
	} elsif ($arg eq "-path") {
	    $showpath = 1;
	} elsif ($arg eq "-cts") {
	    $cts = 1;
	} elsif ($arg eq "-id") {
	    $idcount = 1;
	} else {
	    die $usage, "Unknown option: $arg\n";
	}
    } else {
	push @argv, $arg;
    }
}
die $usage if @argv;

die "-id can only be used with -cts\n" if $idcount && !$cts;

for (my $delrate = $dmin; $delrate <= $dmax; $delrate += $dstep) {

    my $delopen = 1 - exp (-$delrate);

    for (my $gapsize = $gmin; $gapsize <= $gmax; $gapsize += $gstep) {

	my $delext = 1 - 1 / (1 + $gapsize);
	my $insext = $kappa * $delext;

	# insrate = sum_{n=1}^infty (delrate * delext^{n-1} * (1 - delext) * kappa^n) = delrate * kappa * (1 - delext) / (1 - insext)
	my $insrate = $delrate * $kappa * (1 - $delext) / (1 - $insext);
	my $insopen = 1 - exp (-$insrate);

	my %count;  # track transitions between seq lengths
	for (my $rep = 0; $rep < $reps; ++$rep) {

	    warn "[d=$delrate g=$gapsize rep=$rep]\n" if $rep % 1000 == 0;

	    my (@len, @col);
	    my $path = "";
	    if ($cts) {

		# generate ancestral sequence
		my $origlen = $anclen;
		if (!defined $origlen) {
		    $origlen = 0;
		    while (randprob() < $kappa) { ++$origlen }
		}

		my @seq = 0..$origlen-1;
		my $seqlen = $origlen;

		# simulate accumulation of indels as a continous-time process
		my $t = 1;
		while ($t > 0) {
#		    warn "t=$t: (", join (" ", map(defined($_) ? $_ : "i", @seq)), ")";

		    # Let d_k = rate of deletion at position n = L-k  (where 0<=n<L and k = # of sites that can be deleted at position n)
		    # d_k = sum_{n=1}^k x^{k-1} (1-x) = (1-x) \sum_{n=0}^{k-1} x^k = (1-x) S_k
		    #     where x S_k = S_k - 1 + x^k
		    # therefore S_k = (1 - x^k) / (1 - x)
		    #        so d_k = 1 - x^k
		    #      thus d_n = 1 - x^{L-n}
		    my @delratebypos = map ($delrate * (1 - $delext ** ($seqlen-$_)), 0..$seqlen-1);
		    my $totaldel = sum (@delratebypos);
		    my $totalins = $insrate * ($seqlen + 1);
		    my $totalrate = $totaldel + $totalins;
#		    warn "totalins=$totalins delratebypos=(@delratebypos)";

		    $t -= randtime ($totalrate);
		    if ($t > 0) {
			if (randprob() < $totaldel / $totalrate) {

			    # deletion
			    my $delpos = randsample (@delratebypos);
			    my @delratebylen = map ($delext ** $_, 1..$seqlen-$delpos);
			    my $dlen = 1 + randsample (@delratebylen);

			    ++$count{$seqlen}->{$seqlen-$dlen};

			    splice (@seq, $delpos, $dlen);
			    $seqlen -= $dlen;

#			warn "del $dlen at $delpos";

			} else {
			    # insertion
			    my $ilen = 1;
			    while (randprob() < $insext) { ++$ilen }
			    my $inspos = int (rand ($seqlen + 1));

			    ++$count{$seqlen}->{$seqlen+$ilen};

			    splice (@seq, $inspos, 0, map (undef, 1..$ilen));
			    $seqlen += $ilen;

#			warn "ins $ilen at $inspos";

			}
		    }
		}

		# deduce the state path
		my $lastancpos = -1;
		for my $descpos (0..$seqlen-1) {
		    my $ancpos = $seq[$descpos];
		    if (defined $ancpos) {
			while ($ancpos > $lastancpos + 1) {
			    push @col, [++$lastancpos, $descpos];
			    $path .= "d";
			}
			push @col, [$ancpos, $descpos];
			$path .= "m";
			$lastancpos = $ancpos;
		    } else {
			push @col, [$lastancpos, $descpos];
			$path .= "i";
		    }
		}
		while ($origlen > $lastancpos + 1) {
		    push @col, [++$lastancpos, $seqlen-1];
		    $path .= "d";
		}

		@len = ($origlen, $seqlen);

	    } else {

		# simulate state path directly
		@len = (0, 0);
		@col = ([@len]);
		my $deleting = 0;
		while (1) {
		    last if defined($anclen) && $len[0] >= $anclen;
		    last if !defined($anclen) && randprob() >= $kappa;
		    ++$len[0];
		    $deleting = 1 if !$deleting && randprob() < $delopen;
		    if ($deleting) {
			$path .= "d";
		    } else {
			$path .= "m";
			++$len[1];
		    }

		    push @col, [@len];

		    if (!$deleting && randprob() < $insopen) {
			do {
			    $path .= "i";
			    ++$len[1];
			    push @col, [@len];
			} while (randprob() < $insext);
		    }

		    # the following line comes *after* the insertion, *preventing* D->I transitions
		    # (otherwise, deviation peaks at around $delopen ~ 0.5, after which all D's are effectively followed by compensatory I's)
		    $deleting = 0 if $deleting && randprob() >= $delext;
		}
	    }

	    # run-length encode the path
	    for my $state (qw(i m d)) { $path =~ s/($state{3,})/$state@{[length($1)]}/g }

	    # find max deviations
	    my $maxdev = 0;  # max deviation from corner-to-corner diagonal
	    my $maxskew = 0;  # max deviation from y=x
	    if ($len[0] == 0) {   # if $len[0]==0, then can't compute $diag, so compute $maxskew and $maxdev directly
		$maxskew = $maxdev = $len[1];
	    } else {
		for my $pos (@col) {
		    # find intersection of diagonal with 1-coordinate
		    my $diag = $$pos[0] * $len[1] / $len[0];
		    my $dev = abs ($$pos[1] - $diag);
		    my $skew = abs ($$pos[1] - $$pos[0]);

#	warn "delopen=$delopen insopen=$insopen pos=(@$pos) len=(@len) diag=$diag dev=$dev";

		    $maxdev = max ($dev, $maxdev);
		    $maxskew = max ($skew, $maxskew);
		}
	    }

	    my $anc_minus_desc = $len[1] - $len[0];
	    my ($lmin, $lmax) = sort { $a <=> $b } @len;

	    $maxdev = int ($maxdev + .5);  # round up
	    my @path = $showpath ? ($path) : ();
	    print join ("\t", $delrate, $gapsize, $delopen, $insopen, $delext, $insext, $lmax, $anc_minus_desc, $maxdev, $maxskew, @path), "\n";
	}

	# print count of ins->del
	if ($idcount) {
	    for my $n (sort {$a<=>$b} keys %count) {
		for my $m (sort {$a<=>$b} keys %{$count{$n}}) {
		    my $n2m = $count{$n}->{$m};
		    if (exists($count{$m}) && exists($count{$m}->{$n})) {
			my $m2n = $count{$m}->{$n};
			print "IDCOUNT: count($n->$m) = $n2m, count($m->$n) = $m2n, ratio = ", $n2m/$m2n, ", expected ", $kappa ** ($m - $n), "\n";
		    } else {
			print "IDCOUNT: count($n->$m) = $n2m\n";
		    }
		}
	    }
	}
    }
}

# max
sub max {
    my ($max, @x) = @_;
    for my $x (@x) { $max = $x if $x > $max }
    return $max;
}

# random number from 0 to 1
sub randprob {
    return rand(1);
}

# random wait time
sub randtime {
    my ($rate) = @_;
    return -log(randprob()) / $rate;
}

# sample from array of (unnormalized) probabilities
sub randsample {
    my @prob = @_;
    return undef unless @prob;
    my ($rand, $n);
    for ($n = 0, $rand = rand (sum (@prob)); $rand > $prob[$n]; $rand -= $prob[$n++]) { }
    return $n;
}

# sum a series
sub sum {
    my @x = @_;
    my $total = 0;
    grep ($total += $_, @x);
    return $total;
}
