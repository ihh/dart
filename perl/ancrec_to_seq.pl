#!/usr/bin/perl -w

use Getopt::Long;
use Stockholm;

my $tag = "ancrec_CYK_MAP";
my $keepgs = 0;
my $usage = "$0 [-tag <TAG>] [-keepgs] <file>\nCopies '#=GR TAG SEQID ...' lines to 'SEQID ...' lines\nDefault TAG is $tag\n";

GetOptions ("tag=s" => \$tag,
	    "keepgs" => \$keepgs) or die $usage;

@ARGV = qw(-) unless @ARGV;
for my $filename (@ARGV) {
    my $stock = Stockholm->from_file ($filename);
    for my $seqname (@{$stock->seqname}) {
	my @gr = $stock->gr->{$tag}->{$seqname};
	if (@gr > 0) {
	    $stock->seqdata->{$seqname} = join ("", @gr);
	    delete $stock->gr->{$tag}->{$seqname};
	}
    }
    unless ($keepgs) { $stock->{'gs'} = {} }
    print $stock->to_string;
}

