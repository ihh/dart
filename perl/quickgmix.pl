#!/usr/bin/perl -w

my $usage = "\nUsage: $0 <number of mixture components> [file]\n\nUses EM to estimate ML fit of mixture-geometric to empirical nonzero-length distribution\nExpects one- or two-column input format: LENGTH [FREQUENCY]\n";

die $usage if @ARGV < 1 || @ARGV > 2 || grep /^-/, @ARGV;

# init model
my $cpts = shift @ARGV;
die "$usage\nNeed an integer number of components\n" unless $cpts == int($cpts) && $cpts >= 1;

my $seed_base = 0.5;  # range of initial param values is evenly spaced between $seed_base and 1.0
my $pseud_avlen = 0.1;  # pseudocount for avlen statistic (see below)
my $pseud_total = 0.1;  # pseudocount for total statistic (see below)

my @weight = map (1/$cpts, 1..$cpts);
my @param = map ($seed_base + (1-$seed_base)*$_/($cpts+1), 1..$cpts);

# read data
my %freq;
my $line = 1;
while (<>) {
    if (/\S/) {
	my @f = split;
	push @f, 1 if @f == 1;
	die "Bad input format at line $line -- should be one column (LENGTH) or two columns (LENGTH FREQUENCY)\n" if @f != 2;
	die "Bad input format at line $line -- first column (LENGTH) should be a positive integer\n" if $f[0] < 1 || $f[0] != int($f[0]);
	$freq{$f[0]} += $f[1];
    }
    ++$line;
}

# do EM
my $prev_loglike;
for (my $iter = 1; ; ++$iter) {

    # Suppose y[n] is the length of sample n, then the model is...

    #   P(y[n]) = sum_i P(y[n],i)
    # P(y[n],i) = weight[i] * param[n]^(y[n]-1) * (1-param[n])

    # EM equations lead to the following....

    #   param[n] <-- (sum_i c_i[n] (y[n] - 1)) / (sum_i c_i[n] y[n])
    #  weight[n] <-- (sum_i c_i[n]) / (sum_k sum_i c_i[k])

    # ...where c_i[n] is the posterior probability that y[n] came from component i...

    #    c_i[n] = P(y[n],i) / P(y[n])

    # If y[n] is observed f[n] times, these equations become...

    #   param[n] <-- (sum_i c_i[n] f[n] (y[n] - 1)) / (sum_i c_i[n] f[n] y[n])
    #  weight[n] <-- (sum_i f[n] c_i[n]) / (sum_k sum_i f[k] c_i[k])

    # So the summary statistics are...

    # avlen[n] = sum_i c_i[n] f[n] y[n]
    # total[n] = sum_i c_i[n] f[n]

    # ...and in terms of these summary statistics...

    #  param[n] <-- (avlen[n] - total[n]) / avlen[n]
    # weight[n] <-- total[n] / (sum_k total[k])

    # calculate summary statistics and total log-likelihood (E-step)
    my $loglike = 0;
    my @avlen = map ($pseud_avlen, 1..$cpts);
    my @total = map ($pseud_total, 1..$cpts);
    while (my ($len, $freq) = each %freq) {
	# calculate posterior probabilities
	my @cpt_loglike = map (log($weight[$_]) + log($param[$_]) * ($len - 1) + log(1-$param[$_]), 0..$cpts-1);
	my $max_cpt_loglike = max (@cpt_loglike);

	$loglike += $freq * $max_cpt_loglike;
	@cpt_loglike = map ($_ - $max_cpt_loglike, @cpt_loglike);

	my $prob = 0;
	for my $cpt_loglike (@cpt_loglike) {
	    $prob += exp($cpt_loglike);
	}
	$loglike += $freq * log($prob);

	my @postprob = map (exp($cpt_loglike[$_]) / $prob, 0..$cpts-1);
	@avlen = map ($avlen[$_] + $postprob[$_] * $freq * $len, 0..$cpts-1);
	@total = map ($total[$_] + $postprob[$_] * $freq, 0..$cpts-1);

	# debug message
#	warn "loglike=$max_cpt_loglike+(@cpt_loglike) postprob=(@postprob)";
    }

    # print log message
    warn "Iteration $iter: log-likelihood = $loglike, weight=(@weight), param=(@param), avlen=(@avlen), total=(@total)\n";
    last if defined($prev_loglike) && $loglike <= $prev_loglike;
    $prev_loglike = $loglike;

    # update (M-step)
    my $weight_norm = 0;
    for my $total (@total) { $weight_norm += $total }
    @param = map (($avlen[$_] - $total[$_]) / $avlen[$_], 0..$cpts-1);
    @weight = map ($total[$_] / $weight_norm, 0..$cpts-1);

    # sort by param
    my @schwartz = sort { $param[$a] <=> $param[$b] } 0..$cpts-1;
    @param = @param[@schwartz];
    @weight = @weight[@schwartz];
}

# output final params
for my $cpt (0..$cpts-1) {
    print "Component ", $cpt+1, "\tWeight $weight[$cpt]\tParam $param[$cpt]\n";
}


# maximum of a list
sub max {
    my ($max, @x) = @_;
    for my $x (@x) {
	$max = $x if $x > $max;
    }
    return $max;
}
