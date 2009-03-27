#!/usr/bin/perl -w

use strict;

my @tests = qw /simpleTest bifurcTest tkfst/;

my $testfile = "999.tm";
foreach my $test (@tests) {
  my $singletfile = "../tt/$test.singlet.tt";
  my $branchfile = "../tt/$test.branch.tt";

  # valid TM files to compare against
  my $tmfile = "../tt/$test.valid.tm";
  my $rtmfile = "../tt/$test.valid.reduced.tm";
  my $frtmfile = "../tt/$test.valid.fullyreduced.tm";

  # TM
  warn "\n**** testing $test ****";
  warn "checking TM...\n";
  system "./fourWay.pl --tm $singletfile $branchfile > $testfile";
  system "./compareTM.pl $tmfile $testfile";
  unless ($? eq 0) { die "Test failed on $test.\n"; } # $? gets the exit status of the last system call

  # reduced TM
  warn "checking reduced TM...\n";
  system "./fourWay.pl --rtm $singletfile $branchfile > $testfile";
  system "./compareTM.pl $rtmfile $testfile";
  unless ($? eq 0) { die "Test failed on $test.\n"; }

  # fully-reduced TM
  warn "checking fully-reduced TM...\n";
  system "./fourWay.pl --frtm $singletfile $branchfile > $testfile";
  system "./compareTM.pl $frtmfile $testfile";
  unless ($? eq 0) { die "Test failed on $test.\n"; }
}

system "rm $testfile";

warn "Tests out ok.\n\n";

exit 0;
