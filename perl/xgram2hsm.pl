#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use PhyloGram;

# Convert xgram-format files into old-school HSM files (input to tkfemit)

# % strings dart/bin/xrate       (or cat dart/src/*/*)
#  | grep -A20 "File format for hidden substitution matrix"
# File format for hidden substitution matrix
# ==========================================
# Line 0 has two integer parameters "C A"
#  where C is the number of hidden classes
#        A is the number of alphabet symbols (e.g. 20 for proteins)
# Line 1 has A whitespace-separated alphabet tokens
# Line 2 has C*A floating-point parameters representing probabilities at equilibrium
# Lines    3..A+2     are the intra-class substitution matrix for class #1
#                      (A*A reversible rate matrix, R_ij = rate from i to j)
# Lines  A+3..2A+2    are the intra-class matrix for class #2
#   ...etc...
# Lines    X..X+C-1  (where X=CA+3) are the inter-class matrix for residue #1
# Lines  X+C..X+2C-1   are the inter-class matrix for residue #2
# All numeric parameters are delimited by whitespace.

my $xgramFile = shift || "-";
my $gram = PhyloGram->from_file ($xgramFile);

my @token = @{$gram->alphabet->token->value};

my @chains = $gram->all_chains;
die "More than one chain" if @chains > 1;

my $chain = $chains[0];
my @term = @{$chain->terminal->value};
die "More than one terminal" if @term > 1;

my @class = $chain->find_all("hidden-class") ? @{$chain->hidden_class->label->value} : (1);

print @class+0, " ", @token+0, "\n";
print "@token\n";

my @p;
for my $c (0..@class-1) {
    for my $i (0..@token-1) {
	my $p = $chain->initial (state($i,$c));
	push @p, defined($p) ? $p->prob->value : 0.;
    }
}
print "@p\n";

my $mutateHash = $chain->mutate_hash;
#while (($k,$v) = each %$mutateHash) { print "$k $v\n" }

for my $c (0..@class-1) {
    for my $i (0..@token-1) {
	my @r = map (0, @token);
	for my $j (0..@token-1) {
	    if ($j != $i) {
		my $r = $chain->mutate (state($i,$c), state($j,$c), $mutateHash);
		if (defined $r) {
		    $r[$j] = $r->rate->value;
		    $r[$i] -= $r[$j];
		}
	    }
	}
	print "@r\n";
    }
}

for my $i (0..@token-1) {
    for my $c (0..@class-1) {
	my @r = map (0, @class);
	for my $d (0..@class-1) {
	    if ($d != $c) {
		my $r = $chain->mutate (state($i,$c), state($i,$d), $mutateHash);
		if (defined $r) {
		    $r[$d] = $r->rate->value;
		    $r[$c] -= $r[$d];
		}
	    }
	}
	print "@r\n";
    }
}


sub state {
    my ($i, $c) = @_;
    return @class>1 ? [$token[$i],$class[$c]] : [$token[$i]];
}
