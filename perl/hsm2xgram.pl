#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

# Convert old-school HSM files (input to tkfemit) into xgram-format files
my $alphName = "ALPHABET";
my $gramName = "GRAMMAR";
my $chainName = "MATRIX";
my $termName = "EMISSION";
my $stateName = "STATE";
my $className = "CLASS";
my $emitSuffix = "*";

# Read number of classes & alphabet size
$_ = <>;
my ($classes, $alphSize) = split;
my $states = $classes * $alphSize;

# Read alphabet
$_ = <>;
my @alph = split;
die "Alphabet size is ", @alph+0, "; expected $alphSize" if $alphSize != @alph;

# Read initial distribution
$_ = <>;
my @init = split;
die "Initial distribution line has ", @init+0, " entries; expected $states" if $states != @init;

# Read intra- and inter-class matrices
my @intra;
for (my $c = 0; $c < $classes; ++$c) {
    my $intra = {};
    for (my $i = 0; $i < $alphSize; ++$i) {
	$_ = <>;
	my @row = split;
	die "Row $i of intra-class matrix \#$c has ", @row, " entries; expected $alphSize" if $alphSize != @row;
	$intra->{$i} = {};
	for (my $j = 0; $j < $alphSize; ++$j) {
	    $intra->{$alph[$i]}->{$alph[$j]} = $row[$j] if $i != $j;
	}
    }
    push @intra, $intra;
}

my %inter;
for my $i (@alph) {
    my $inter = {};
    for (my $c = 0; $c < $classes; ++$c) {
	$_ = <>;
	my @row = split;
	die "Row $c of inter-class matrix for residue $alph[$i] has ", @row, " entries; expected $classes" if $classes != @row;
	for (my $d = 0; $d < $classes; ++$d) {
	    $inter->{$c}->{$d} = $row[$d] if $c != $d;
	}
    }
    $inter{$i} = $inter;
}

# Read anything at the end
my %tagval;
while (<>) {
    my ($tag, $val) = split;
    $tagval{$tag} = $val;
}

# Set some default tag-values
sub set_tagval { my ($tag, $val) = @_; $tagval{$tag} = $val unless exists $tagval{$tag} }
set_tagval ('update-policy', 'rev');
set_tagval ('update-rates', 1);
set_tagval ('update-rules', 0);

# Display
print ";; The alphabet\n";
print "(alphabet (name $alphName) (token (@alph)) (wildcard *))\n\n";

print ";; The grammar\n";
print "(grammar\n";
print " (name $gramName)\n";
print " (update-rates ", $tagval{'update-rates'}, ")\n";
print " (update-rules ", $tagval{'update-rules'}, ")\n";

# Print grammar rules
print " ;; Null model grammar rules: emit, self-loop, terminate\n";
print " (transform (from ($stateName)) (to ($termName $stateName$emitSuffix)))\n";
print " (transform (from ($stateName$emitSuffix)) (to ($stateName)))\n";
print " (transform (from ($stateName$emitSuffix)) (to ()))\n\n";

# Print chain
print " ;; Rate matrix\n";
print " (chain\n";
print "  (update-policy ", $tagval{'update-policy'}, ")\n";
print "  (terminal ($termName))\n";
print "  (hidden-class (row $className) (label (", join(' ',0..$classes-1), ")))\n" if $classes > 1;

print "\n";
print "  ;; Initial distribution over states\n";
for (my $c = 0; $c < $classes; ++$c) {
    for (my $i = 0; $i < @alph; ++$i) {
	print "  (initial (state (", state($alph[$i],$c,$classes), ")) (prob ", $init[$c * $alphSize + $i], "))\n";
    }
}

print "\n";
print "  ;; Mutation rates\n";
for (my $c = 0; $c < $classes; ++$c) {
    for my $i (@alph) {
	for my $j (@alph) {
	    print "  (mutate (from (", state($i,$c,$classes), ")) (to (", state($j,$c,$classes), ")) (rate ", $intra[$c]->{$i}->{$j}, "))\n" if $j ne $i;
	}
	for (my $d = 0; $d < $classes; ++$d) {
	    print "  (mutate (from (", state($i,$c,$classes), ")) (to (", state($i,$d,$classes), ")) (rate ", $inter{$i}->{$c}->{$d}, "))\n" if $d != $c;
	}
    }
}
print " )  ;; end of rate matrix\n";
print ")  ;; end of grammar\n";


# Subroutine to print a state
sub state {
    my ($i, $c, $classes) = @_;
    return $classes > 1 ? "$i $c" : $i;
}
