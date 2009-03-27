#!/usr/bin/perl

my @x;
while (<STDIN>) { push @x, $_ }

warn "\n";
warn "@ARGV" if @ARGV;
warn "Start of standard input";
warn @x if @x;
warn "End of standard input";
warn "\n";

print STDOUT map ("$_ ", @ARGV);
print STDOUT @x;
