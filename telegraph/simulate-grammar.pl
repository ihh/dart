#!/usr/bin/perl -w

# init
my $DELAY = 10**1;
my $ENDSTATE = "end";
my $VAR = "[A-Za-z_][A-Za-z_0-9]*";
my $LHS = "$VAR( $VAR)*";
my $RHS = "|$VAR( $VAR)*";
my $STRING = "'[^']+'|\"[^\"]+\"";
my $NUM = '[0-9\.\+\-\*\/\(\)eE]+';
my $RULERATE = "\{ ?$NUM ?\}";

my $usage = "Grammar simulator (Telegraph subset)\n";
$usage .= "Usage: $0 [<file>]\n";

# parse args
my @argv;
while (@ARGV) {
    my $argv = shift @ARGV;
    if ($argv =~ /^-/) {
	die "Unrecognised option $argv\n\n", $usage;
    } else {
	push @argv, $argv;
    }
}

# read input
die $usage unless @argv < 2;
@ARGV = @argv;
$/ = ";";  # split on semicolons
my @line = <>;
grep (chomp, @line);  # remove semicolons
$/ = "\n";

# remove "end" keywords
grep (s/\b$ENDSTATE\b//g, @line);

# remove unnecessary whitespace
grep (s/^\s*//, @line);
grep (s/\s*$//, @line);
grep (s/\s+/ /g, @line);

# parse grammar
my %rules;
my $maxRhsLen = 0;
my $start;
foreach (@line) {

    next unless /\S/;

    # match "LHS -> ..."
    if (/^(.*?) ?-> ?(.*)$/) {
	my ($lhs, $rhsBlock) = ($1, $2);

	# warn "lhs=$lhs rhsBlock=$rhsBlock";
	die "Syntax error: $_\n" unless $lhs =~ /^$LHS$/;
	$start = $lhs unless defined $start;

	# match "RHS1 { RATE1 } | RHS2 { RATE2 } ..."
	while (1) {
	    if ($rhsBlock =~ /^(.*?) ?(|$RULERATE) ?\| ?(.*)$/) {  # multiple RHS's
		my ($rhs, $ruleRate, $nextRhsBlock) = ($1, $2, $3);
		die "Syntax error: $_\n" unless $rhs =~ /^$RHS$/;
		addRule (\%rules, \$maxRhsLen, $lhs, $rhs, $ruleRate);
		$rhsBlock = $nextRhsBlock;

	    } elsif ($rhsBlock =~ /^(.*?) ?(|$RULERATE) ?$/) {  # single RHS
		my ($rhs, $ruleRate) = ($1, $2);
		die "Syntax error: $_\n" unless $rhs =~ /^$RHS$/;
		addRule (\%rules, \$maxRhsLen, $lhs, $rhs, $ruleRate);
		last;

	    } else {  # syntax error
		die "Syntax error: $_\n";
	    }
	}
    }
}

# set output autoflush to true
$| = 1;

# check we've got a start nonterm
die "Need at least one nonterminal\n" unless defined $start;

# simulate
my @seq = ($start);
while (1) {
    # print current sequence
    my $newline = 1;
    foreach (@seq) {
	print $newline ? "\n" : " ", $_;
	$newline = 0;

	# delay loop
	for ($i = 0; $i < $DELAY; ++$i) {}
    }

    # find all applicable rules
    my $totalRate = 0;
    my @rule;
    for (my $len = 1; $len <= $maxRhsLen; ++$len) {
	for (my $pos = 0; $pos <= @seq - $len; ++$pos) {
	    my $lhs = join (" ", @seq[$pos..$pos+$len-1]);
	    if (exists $rules{$lhs}) {
		foreach my $rhs_rate (@{$rules{$lhs}}) {
		    my ($rhs, $rate) = @$rhs_rate;
		    $totalRate += $rate;
		    push @rule, [$rate, $pos, $len, \$rhs];
		    # warn "[applicable rule at $pos+$len: $lhs -> $rhs { $rate }]\n";
		}
	    }
	}
    }

    # exit if terminated
    last unless $totalRate;

    # sample a rule and apply it
    my $r = rand() * $totalRate;
    foreach my $rule (@rule) {
	my ($rate, $pos, $len, $rhs) = @$rule;
	if (($r -= $rate) <= 0) {
	    # warn "[applying rule at $pos+$len: ", join(" ",@seq[$pos..$pos+$len-1]), " -> $$rhs]\n";
	    splice (@seq, $pos, $len, split (/ /, $$rhs));
	    last;
	}
    }
}

# addRule method
sub addRule {
    my ($rulesHashRef, $maxRhsLenRef, $lhs, $rhs, $rate) = @_;
    if ($rate eq "") {
	$rate = 1;
    } elsif ($rate =~ /($NUM)/) {
	$rate = $1;
    } else {
	die "Bad rate: $rate\n";
    }

    # warn "Adding rule: $lhs -> $rhs { $rate };\n";
    if (!exists $rulesHashRef->{$lhs}) {
	$rulesHashRef->{$lhs} = [];
    }
    push @{$rulesHashRef->{$lhs}}, [$rhs, $rate];
    my @rhs = split / /, $rhs;
    $$maxRhsLenRef = @rhs if @rhs > $$maxRhsLenRef;
}
