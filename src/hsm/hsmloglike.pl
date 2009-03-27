#!/usr/local/bin/perl -w

push @ARGV, '-' unless @ARGV;
foreach my $file (@ARGV) {
    print "\"", ($file eq '-' ? "STDIN" : $file), "\n";
    local *FILE;
    open FILE, "<$file" or die $!;
    while (<FILE>) {
	if (/EM iteration \#(\d+): log-likelihood = (\S+)/) {
	    my ($iter, $loglike) = ($1, $2);
	    print "$iter $loglike\n";
	}
    }
    close FILE;
    print "\n";
}
