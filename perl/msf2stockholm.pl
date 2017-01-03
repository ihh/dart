#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

my $width = 50;

my $found_sep = 0;
while (<>) { if (/^\/\//) { $found_sep = 1; last } }
die "Couldn't find // separator" unless $found_sep;

my @seqname;
my %seq;
while (<>) {
    next if /^\s*\#/ || !/\S/;
    my @f = split /\s+/, $_, 2;
    $f[1] =~ s/\s//g;
    $f[1] =~ s/-/./g;
    if (!exists $seq{$f[0]}) { push @seqname, $f[0] }
    $seq{$f[0]} .= $f[1];
}

print "# STOCKHOLM 1.0\n";
my $nwidth = max (map (length($_), keys %seq));
while (my ($name,$seq) = each %seq) {
    if ($name =~ /\S/) {
	printf "%-".$nwidth."s %s\n", $name, $seq;
    }
}
print "//\n";

sub max {
    my ($max, @x) = @_;
    for my $x (@x) {
	$max = $x if $x > $max;
    }
    return $max;
}
