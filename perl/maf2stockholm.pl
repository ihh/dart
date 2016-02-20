#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm;

# params
my $cols = 50;

# help message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname [-c <column width>] <MAF file>\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-c") {
	defined ($cols = shift) or die $usage;
    } elsif ($arg =~ /^-/) {
	die $usage, "\nUnrecognized argument: $arg\n";
    } else {
	push @argv, $arg;
    }
}
unless (@argv) {
    @argv = ('-');
    warn "[waiting for alignments on standard input]\n";
}
@ARGV = @argv;

my $stock;
my $hadEmpty = 0;
my $isFirst = 1;
while (<>) {
    next if /^\s*\#/;  # comment
    if (/\S/) {
	if (/^a/) {
	    if (!$hadEmpty && !$isFirst) { warn "A new alignment block beginning with an 'a' line should be preceded by an empty line\n" }
	    $hadEmpty = 0;
	    print $stock->to_string($cols) if defined $stock;
	    $stock = Stockholm->new;
	    while (/\s(\S+)=(\S+)/g) { $stock->add_gf ($1, $2) }
	} else {
	    if ($hadEmpty) { warn "Following blank line, expected new alignment block beginning with an 'a' line\n" }
	    $hadEmpty = 0;
	    if (/^s /) {
		my ($s_char, $name, $start, $len, $strand, $end, $seq) = split;
		$stock->add_row ($name, $seq);
		$stock->gs->{'start'}->{$name} = [$start];
		$stock->gs->{'len'}->{$name} = [$len];
		$stock->gs->{'strand'}->{$name} = [$strand];
		$stock->gs->{'end'}->{$name} = [$end];
	    } elsif (/^q/) {
		my ($q_char, $name, $quality) = split;
		$stock->gr->{'quality'}->{$name} = $quality;
	    } elsif (/^[ie]/) {
		# Silently skip MAF 'i' and 'e' lines
	    } elsif (/^track/ && $isFirst) {
		# Silently skip initial 'track' line
	    } else {
		warn "Skipping line: $_";
	    }
	}
    } else {
	$hadEmpty = 1;
    }
    $isFirst = 0;
}
print $stock->to_string($cols) if defined $stock;
