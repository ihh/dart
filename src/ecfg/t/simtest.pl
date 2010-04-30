#!/usr/bin/perl -w

use PhyloGram;
use Stockholm;

die "Usage: $0 [<tree files>]\n" if grep (/^-/, @ARGV);

my $alignSuffix = ".stock";
my $gramSuffix = ".eg";
my $logSuffix = ".log";
my $treeSuffix = ".nh";

my $seedSuffix = ".seed";
my $metaSuffix = ".meta";
my $countsSuffix = ".counts";
my $trainedSuffix = ".trained";

my $perldir = "../../../perl";
my $bindir = "../../../bin";

my $xrate = "$bindir/xrate";
my $simgram = "$bindir/simgram";
my $compareParams = "$perldir/compareParams.pl";

my $logopts = "-log 6";

local *DIR;
opendir DIR, '.';
my @dir = grep -f $_, readdir DIR;
closedir DIR;

my @prefix;
my %trees;
for my $dir (@dir) {
    if ($dir =~ /^(\S+)$gramSuffix$/) {
	my $prefix = $1;
	my @trees = map (/^($prefix.*)$treeSuffix$/ ? ($1) : (), @dir);
	if (@trees) {
	    push @prefix, $prefix;
	    $trees{$prefix} = \@trees;
	}
    }
}

my %gparam_default = qw(initial-tol .1
			mutate-tol .1
			wait-tol .1
			param-tol .1
			max-errs .1
			n-align 1
			min-inc .001
			forgive 0);

my $ok = 1;
foreach my $prefix (@prefix) {
    # get filenames
    my ($gramFile, $baseMetaFile) = map ("$prefix$_", $gramSuffix, "$metaSuffix$gramSuffix");
    warn "Considering grammar $gramFile\n";

    for my $tree (@{$trees{$prefix}}) {
	# get filenames
	my ($treeFile, $alignFile, $metaFile, $treeGramFile, $seedFile, $countsFile, $trainedFile, $logFile) = map ("$tree$_", $treeSuffix, $alignSuffix, "$metaSuffix$gramSuffix", $gramSuffix, "$seedSuffix$gramSuffix", "$countsSuffix$gramSuffix", "$trainedSuffix$gramSuffix", "$trainedSuffix$logSuffix");

	# check @ARGV for explicit treefiles
	next if @ARGV && !grep ($_ eq $treeFile, @ARGV);

	# load grammar & paste metafile(s) into grammar, if it exists
	my $gram = load_meta ($gramFile, $baseMetaFile, $metaFile);
	$gram->to_file ($treeGramFile);

	# set the parameters to the values given in the "seed" block, and save to the ".seed" grammar
	save_seed ($gram, $seedFile);

	# get the tolerances for counts, wait times & trained parameters, and the max error fraction
	my %gparam = get_gparam (\%gparam_default, $gram);

	# reset error count
	my $errs = 0;
	my $tests = 0;

	# simulate an alignment
	syswarn ("$simgram $logopts -g $treeGramFile -t $treeFile >$alignFile");

	my $align = Stockholm->from_file ($alignFile);
	my $real_counts = DartSexpr->from_string (@{$align->gf_("observed-chain-counts")});

	# run xrate using the original grammar for one round of EM, and compare the chain counts
	syswarn ("$xrate $logopts $alignFile -g $treeGramFile -t $countsFile --maxrounds 1 --noannotate");

	my $counts_gram = PhyloGram->from_file ($countsFile);
	my $xrate_counts = $counts_gram->grammar->observed_chain_counts;

	for my $real_counts_block (@{$real_counts}) {
	    my $real_term = "@{$real_counts_block->terminal->values}";

	    my @real_initial = $real_counts_block->find_all ('initial');
	    my @real_mutate = $real_counts_block->find_all ('mutate');
	    my @real_wait = $real_counts_block->find_all ('wait');

	    for my $xrate_counts_block (@{$xrate_counts}[1..@$xrate_counts-1]) {
		my $xrate_term = "@{$xrate_counts_block->terminal->values}";
		if ($xrate_term eq $real_term) {

		    my @xrate_initial = $xrate_counts_block->find_all ('initial');
		    my @xrate_mutate = $xrate_counts_block->find_all ('mutate');
		    my @xrate_wait = $xrate_counts_block->find_all ('wait');

		    my %xrate_mutate_hash;
		    for my $sexpr (@xrate_mutate) {
			$xrate_mutate_hash{join ("",@{$sexpr->from->values},@{$sexpr->to->values})} = $sexpr;
		    }

		    for my $real_initial (@real_initial) {
			my ($xrate_initial) = grep (join("",@{$_->state->values}) eq join("",@{$real_initial->state->values}), @xrate_initial);
			my $cmp = compare ($real_initial->count->value, $xrate_initial->count->value, $gparam{'initial-tol'}, "terminal ($real_term) state (@{$real_initial->state->values}) initial count");
			++$tests;
			++$errs unless $cmp;
		    }

		    for my $real_mutate (@real_mutate) {
			my $xrate_mutate = $xrate_mutate_hash{join("",@{$real_mutate->from->values},@{$real_mutate->to->values})};
			my $cmp = compare ($real_mutate->count->value, $xrate_mutate->count->value, $gparam{'mutate-tol'}, "terminal ($real_term) from (@{$real_mutate->from->values}) to (@{$real_mutate->to->values}) mutate counts");
			++$tests;
			++$errs unless $cmp;
		    }

		    for my $real_wait (@real_wait) {
			my ($xrate_wait) = grep (join("",@{$_->state->values}) eq join("",@{$real_wait->state->values}), @xrate_wait);
			my $cmp = compare ($real_wait->time->value, $xrate_wait->time->value, $gparam{'wait-tol'}, "terminal ($real_term) state (@{$real_wait->state->values}) wait times");
			++$tests;
			++$errs unless $cmp;
		    }
		}
	    }
	}

	# simulate N alignments
	my $nalign = $gparam{'n-align'};
	my $nAlignFile = $alignFile;
	if ($nalign > 1) {
	    $nAlignFile = "$tree.x$nalign$alignSuffix";
	    syswarn ("$simgram $logopts -n $nalign -g $treeGramFile -t $treeFile >$nAlignFile");
	}

	# run xrate from the flat seed and compare the trained parameters to the true ones
	unlink $logFile if -e $logFile;
	my $mininc = $gparam{'min-inc'};
	my $forgive = $gparam{'forgive'};
	syswarn ("$xrate $logopts -logcopy $logFile $nAlignFile -g $seedFile -t $trainedFile --mininc $mininc --forgive $forgive --noannotate");
	my $comp_sexpr = DartSexpr->from_string (syswarn ("$compareParams $treeGramFile $trainedFile"));

	for my $pcomp (@$comp_sexpr) {
	    my $pname = $pcomp->param->value;
	    my $cmp = compare (map ($_->def->value, @{$pcomp}[1..@$pcomp-1]), $gparam{'param-tol'}, "param $pname values");
	    ++$tests;
	    ++$errs unless $cmp;
	}

	# check errors
	print STDERR "\n*** TEST RESULTS ***\nFailed $errs of $tests comparisons for $treeFile ... ";
	my $max_errs = $gparam{'max-errs'};
	if ($errs / $tests > $max_errs) {
	    $ok = 0;
	    print STDERR "outside acceptable fraction ($max_errs), NOT OK\n\n";
	} else {
	    print STDERR "within acceptable fraction ($max_errs), ok\n\n";
	}
    }
}

print $ok ? "ok\n" : "not ok\n";

sub syswarn {
    my ($command) = @_;
    warn "\n***** Running '$command' *****\n";
    return `$command`;
}

sub max {
    my $max = shift;
    for my $x (@_) { $max = $x if $x > $max }
    return $max;
}

sub compare {
    my ($a, $b, $tol, $desc) = @_;
    my $abmax = max($a,$b);
    my ($type,$effmax) = $abmax < 1 ? ("absolute",1) : ("relative",$abmax);
    my $error = abs($a-$b) / $effmax;
    my $ok = $error <= $tol;
    warn "Comparing $desc:   ( $a , $b )   ", $ok ? "passed" : "failed (tolerance $tol, $type error $error)", "\n";
    return $ok;
}

sub get_gparam {
    my ($gparam, $gram) = @_;
    my %gparam = %$gparam;
    for my $gp (keys %gparam) {
	if ($gram->grammar->find_all($gp) > 0) {
	    $gparam{$gp} = $gram->grammar->find($gp)->value;
	}
    }
    return %gparam;
}

sub load_meta {
    my ($gramFile, @metaFile) = @_;
    my $gram = PhyloGram->from_file ($gramFile);
    for my $metaFile (@metaFile) {
	if (-e $metaFile) {
	    my $meta = DartSexpr->from_file ($metaFile);
	    for my $metaChild (@$meta) {
		@{$gram->grammar->find_or_add($metaChild->tag)} = @$metaChild;
	    }
	}
    }
    return $gram;
}

sub save_seed {
    my ($gram, $seedFile) = @_;
    # set the parameters to the values given in the "seed" block, and save to the ".seed" grammar
    my $param_hash = $gram->param_hash;
    my $param_defaults = $gram->grammar->seed;
    while (my ($param, $sexpr) = each %$param_hash) {
	my $seed_value = $param_defaults->find($param)->value;
	$sexpr->value ($seed_value);
    }
    local *SEED;
    open SEED, ">$seedFile";
    print SEED $gram->to_string;
    close SEED;
}
