#!/usr/bin/env perl -w

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
    my ($name, $seq, $len, $lastname, $lastlen);
    open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
    while (<FASTA>) {
	if (/^\s*>\s*(\S+)/) {
	    if (defined ($lastname) && $len != $lastlen) {
		# length changed, so start new alignment
		print "//\n# STOCKHOLM 1.0\n";
		%seen = ();  # forget previously-seen sequence names
	    }
	    print "$name $seq\n" if defined $name;
	    $lastname = $name;
	    $lastlen = $len;
	    $name = $1;
	    $len = 0;
	    $seq = "";
	    warn "Duplicate sequence name: $name\n" if $seen{$name};
	    $seen{$name} = 1;
	} else {
	    if (/\S/ && !defined $name) {
		warn "Ignoring: $_";
	    } else {
		s/\s//g;
		$seq .= $_;
		$len += length;
	    }
	}
    }
    close FASTA;
    if (defined ($lastname) && $len != $lastlen) {
	# length changed, so start new alignment
	print "//\n# STOCKHOLM 1.0\n";
    }
    print "$name $seq\n" if defined $name;
    print "//\n";
}
