#!/usr/bin/perl -w
# Pete K begun: 3/15/05
# test_exphsm.pl - invoke ".sym", EM_matrix, and ".asym",
# Irrev_EM_matrix versions of ian's exphsm program,
# which reads in an xrate hsm format file,
# calls create_joint_substitution_matrix, and
# outputs the calculated matrix.
# Tests to run are defined in file test_exphsm.in.
# The first line of test_exphsm.in contains the 
# maximum allowed absolute difference between matrix
# elements (tolerance).  Following lines describe
# tests, one per line - for each the .hsm file
# and the time to use.
# For each test, invoke exphsm, irrev_exphsm,
# writing output to exphsm.tmp and irrev_exphsm.tmp
# Read in exphsm.tmp and irrev_exphsm.tmp arrays,
# calculate the absolute difference between respective array
# elements, print out the maximum difference, and compare
# that with the accepted tolerance (first line of test_exphsm.in).
# If all absolute differences are below the tolerance,
# output 'ok', else 'not ok'.
use strict;

my $test_exphsm_in_file = "test_exphsm.in";
my $exphsm_out_file = "exphsm.tmp";
my $irrev_exphsm_out_file = "irrev_exphsm.tmp";
my $irrevOutSuffix = ".irrev_exphsm.out";

my $exphsm = "./exphsm";
my $irrev_exphsm = "./irrev_exphsm";

my $test_line;
open TESTS_IN, "<$test_exphsm_in_file" or die "'$test_exphsm_in_file': $!";
if (!defined ($test_line = <TESTS_IN>) ) {
    die "no tolerance defined in input file\n";
}
chomp $test_line;
my @test_fields = split /\s+/, $test_line;
if ($#test_fields != 0) {
	die "invalid tolerance line\n";				# expect only one entry, the tolerance
}
my $tol = $test_fields[0];

my $ok = 1;
my $test_count = 1;
while (defined ($test_line = <TESTS_IN>) ) {
	chomp $test_line;
	@test_fields = split /\s+/, $test_line;
	if ($#test_fields != 1) {
		die "test $test_count invalid data line\n";				# two entries per line
	}
	my $hsm_file = $test_fields[0];
	my $time = $test_fields[1];

    my $command = "$exphsm -t $time $hsm_file >$exphsm_out_file";
    warn "[Invoke '$command' ]\n";
    system $command;

    $command = "$irrev_exphsm -t $time $hsm_file >$irrev_exphsm_out_file";
    warn "[Invoke '$command' ]\n";
    system $command;

	my $max_diff = 0;

    open EXPHSM_IN, "<$exphsm_out_file" or die "'$exphsm_out_file': $!";
    open IR_EXPHSM_IN, "<$irrev_exphsm_out_file" or die "'$irrev_exphsm_out_file': $!";
	my ($exphsm_line, $ir_exphsm_line);
	while (defined ($exphsm_line = <EXPHSM_IN>) ) {
		if (!defined ($ir_exphsm_line = <IR_EXPHSM_IN>) ) {
			die "test $test_count line count mismatch between exphsm, irrev_exphsm output\n";
		}
		chomp $exphsm_line;
		chomp $ir_exphsm_line;
		my @exphsm_fields = split /\s+/, $exphsm_line;
		my @ir_exphsm_fields = split /\s+/, $ir_exphsm_line;
		if ($#exphsm_fields != $#ir_exphsm_fields) {
			die "test $test_count field count mismatch between exphsm, irrev_exphsm output\n";
		}
		for (my $i = 0; $i <= $#exphsm_fields; $i++) {
			my $diff = $exphsm_fields[$i] - $ir_exphsm_fields[$i];
			$diff = $diff >= 0 ? $diff : -$diff;
			$max_diff = $diff > $max_diff ? $diff : $max_diff;
		}
	}
    if (defined (my $extra_line = <IR_EXPHSM_IN>) ) {
		die "test $test_count line count mismatch between exphsm, irrev_exphsm output\n";
	}

	my $testok =  $max_diff <= $tol ? 1 : 0;
	$ok = $testok ? $ok : 0;

    warn "test $test_count maximum difference $max_diff\n";
	print "test $test_count ";
	print $testok ? "ok\n" : "not ok\n";
	$test_count++;
}

unlink $exphsm_out_file, $irrev_exphsm_out_file;
print $ok ? "ok\n" : "not ok\n";
