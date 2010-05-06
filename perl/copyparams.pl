#!/usr/bin/env perl -w

use PhyloGram;

my $progname = $0;
$progname =~ s!^.*/!!;

my $usage = "$progname: copy parameters from one grammar to another\n";
$usage .= "(Parameters must already be defined somewhere in destination grammar)\n";
$usage .= "\n";
$usage .= "Usage: $progname [options] <source-grammar> <dest-grammar>\n";
$usage .= "\n";
$usage .= "Options:\n";
$usage .= " -e <regexp>   Only copy parameters whose name matches the given regular expression\n";
$usage .= "\n";

my $regexp;
my @argv;
while (@ARGV) {
    my $argv = shift;
    if ($argv eq "-e") {
	defined ($regexp = shift) or die $usage;
    } elsif ($argv =~ /^-./) {
	die $usage, "\nUnknown option: $argv\n";
    } else {
	push @argv, $argv;
    }
}

die $usage unless @argv == 2;
my ($srcGram, $destGram) = map (load_or_die($_), @argv);
my ($srcHash, $destHash) = map ($_->param_hash, $srcGram, $destGram);

my (@notInDest, @notMatch);
while (my ($param, $srcExpr) = each %$srcHash) {
    if (defined $regexp) {
	unless ($param =~ /$regexp/) {
	    push @notMatch, $param;
	    next;
	}
    }
    unless (exists $destHash->{$param}) {
	push @notInDest, $param;
	next;
    }
    @{$destHash->{$param}} = @$srcExpr;
}

warn "Skipped the following parameters, as they did not match the regular expression: @{[sort @notMatch]}\n" if @notMatch;
warn "Skipped the following parameters, as they were not defined in the destination grammar: @{[sort @notInDest]}\n" if @notInDest;

print $destGram->tidy_string;

sub load_or_die {
    my ($filename) = @_;
    die "File '$filename' not found" unless -e $filename || $filename eq '-';
    my $grammar = PhyloGram->from_file ($filename);
    die "Failed to load grammar '$filename'" unless defined $grammar;
    return $grammar;
}
