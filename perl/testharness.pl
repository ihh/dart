#!/usr/bin/env perl -w

# Usage: $0 [...tests in dart/bin dir...] -fullpath [-cd 'testdir'] [...tests with fully qualified pathnames...]

my $bindir = $0;
$bindir =~ s/\/[^\/]+$//;
$bindir .= "/../bin";

my $progname = $0;
$progname =~ s/^.*\///;

my @tests;
my $fullpath;
while (@ARGV) {
    my $argv = shift;
    if ($argv eq "-fullpath") {
	$fullpath = 1;
    } elsif ($argv eq "-cd") {
	my $newtestdir = shift;
	chdir $newtestdir or warn "[$progname: can't chdir to '$newtestdir'; ignoring '$argv' directive]\n";
    } else {
	$argv = "$bindir/$argv" unless defined $fullpath;
	push @tests, $argv;
    }
}

my @done;
foreach my $test (@tests) {
    my $brief = $test;
    $brief =~ s/^.*\///;
    print "Commencing test: $brief ($test)\n";
    my @out = `$test`;
    if (!grep (/^\s*ok\s*$/, @out) || grep (/not ok/, @out)) { die map ("> $_", "Test '$test' failed, output:\n", @out), "\nnot ok\n" }
    print "Passed test: $brief ($test)\n";
    push @done, $brief;
}
print "Passed all tests: @done\n";
print "ok\n";
