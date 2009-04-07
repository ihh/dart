#!/usr/bin/perl -w

use Stockholm::Database;

# options
my $delay_secs = .02;  # delay in seconds between "frames"
my $cols;  # alignment columns per row

# help message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname <Stockholm file>\n";
$usage .=   "          [-h] print this message\n";
$usage .=   "      [-delay] frame delay in seconds (default is $delay_secs)\n";
$usage .=   "   [-cols <N>] number of columns per row\n";
$usage .=   "\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-cols") {
	defined ($cols = shift) or die $usage;
    } elsif ($arg eq "-delay") {
	defined ($delay_secs = shift) or die $usage;
    } else {
	push @argv, $arg;
    }
}
unless (@argv) {
    @argv = ('-');
    warn "[waiting for alignments on standard input]\n";
}
die $usage unless @argv == 1;
my ($filename) = @argv;

my $db = Stockholm::Database->from_file ($filename);
warn "Animation will take ", $delay_secs * @$db, " seconds (plus display overhead)\n";

warn "Calibrating delay loop...\n";
my $ticks_per_sec = 0;
alarm 1;
$SIG{'ALRM'} = \&alarm_sub;   # hmm, must be a less hacky way of doing this...
while (1) { ++$ticks_per_sec }

sub alarm_sub {
    alarm 0;
    my $delay_ticks = $delay_secs * $ticks_per_sec;

    for my $stock (@$db) {
	system "clear";  # clear the screen
	print $stock->to_string ($cols);
	for (my $ticks = 0; $ticks < $delay_ticks; ++$ticks) { }
    }

    exit;
}
