#!/usr/bin/perl -w

my $usage = "Usage: $0 <gapped FASTA alignment file(s)>\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-/) {
	if ($arg eq "-h") { print $usage; exit }
	else { die $usage }
    } else {
	push @argv, $arg;
    }
}
push @argv, "-" unless @argv;

# loop through FASTA files
foreach my $fasta (@argv) {
# print Stockholm output
    print "# STOCKHOLM 1.0\n";

# read FASTA file
    my %len;
    my @name;
    my $name;
    open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
    while (<FASTA>) {
	if (/^\s*>\s*(\S+)/) {
	    print "\n" if defined $name;
	    $name = $1;
	    die "Duplicate sequence name: $name\n" if exists $len{$name};
	    print "$name ";
	    push @name, $name;
	} else {
	    if (/\S/ && !defined $name) {
		warn "Ignoring: $_";
	    } else {
		s/\s//g;
		print;
		$len{$name} += length;
	    }
	}
    }
    close FASTA;
    print "\n" if defined $name;
    print "//\n";

    warn map (" $_ : $len{$_}\n", @name);

# check all seqs are same length
    my $length;
    foreach my $name (@name) {
	my $l = $len{$name};
	if (defined $length) {
	    if ($length != $l) {
		die "Sequences not all same length:\n", map (" $_ : $len{$_}\n", @name);
	    }
	} else {
	    $length = $l;
	}
    }
}
