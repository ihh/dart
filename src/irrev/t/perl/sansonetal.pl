#!/usr/bin/perl -w

my @p = ([0.706064, 0.0167902, 0.225862, 0.0512831], [0.0113993, 0.907438,
     0.00193142, 0.0792316], [0.0940831, 0.0119462, 
    0.874856, 0.0191143], [0.0641252, 0.221503, 0.0176224, 0.69675]);

my @a = qw(a c g t);

my ($sx, $sy) = ("", "");
for (my $i = 0; $i < 1000; ++$i) {
	my $p = rand(1);
	my $x = int rand(4);
	for ($y = 0; $y < 4; ++$y) {
		last if ($p -= $p[$x]->[$y]) <= 0;
	}
	$sx .= $a[$x];
	$sy .= $a[$y];
}

print "# STOCKHOLM 1.0\n";
print "#=GF NH (seqx:.01,seqy:1);\n";
print "seqx $sx\nseqy $sy\n";
print "//\n";
