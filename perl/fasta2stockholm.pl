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
    my %seen;
    my ($name, $len, $lastname, $lastlen);
    open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
    while (<FASTA>) {
	if (/^\s*>\s*(\S+)/) {
	    print "\n" if defined $name;
	    if (defined ($lastname) && $len != $lastlen) {
		die "Sequence $name (length $len) has different length from sequence $lastname (length $lastlen)\n";
	    }
	    $lastname = $name;
	    $lastlen = $len;
	    $name = $1;
	    $len = 0;
	    warn "Duplicate sequence name: $name\n" if $seen{$name};
	    $seen{$name} = 1;
	    print "$name ";
	} else {
	    if (/\S/ && !defined $name) {
		warn "Ignoring: $_";
	    } else {
		s/\s//g;
		print;
		$len += length;
	    }
	}
    }
    close FASTA;
    print "\n" if defined $name;
    if (defined ($lastname) && $len != $lastlen) {
	die "Sequence $name (length $len) has different length from sequence $lastname (length $lastlen)\n";
    }
    print "//\n";
}
