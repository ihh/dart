#!/usr/bin/perl -w

# convert RNAfold output to Stockholm

my $name = "AnonymousSequence";
my (%seq, %ss, @name);
while (<>) {
    if (/^\s*>\s*(\S+)/) {
	$name = $1;
	die "Duplicate name '$name'" if exists $seq{$name};
	$seq{$name} = $ss{$name} = "";
	push @name, $name;
    } elsif (defined $name) {
	push @name, $name unless @name > 0;  # handle the case of anonymous sequences
	if (/^\s*(\S+)\s*$/) {
	    $seq{$name} .= $1;
	} elsif (/^\s*([\(\)\.]+)\s+\(\s*[0-9\.\-]+\)/) {
	    my $ss = $1;
	    $ss =~ tr/()/<>/;
	    $ss{$name} .= $ss;
	}
    }
}

for my $name (@name) {
    my $w = max (length($name), 4) + 8;

    print "# STOCKHOLM 1.0\n";
    printf "%-${w}s $seq{$name}\n", $name;
    printf "%-${w}s $ss{$name}\n", "#=GR $name SS";
    printf "%-${w}s $ss{$name}\n", "#=GC SS_cons";
    print "//\n";
}

sub max {
    my ($max, @x) = @_;
    for my $x (@x) { $max = $x if $x > $max }
    return $max;
}
