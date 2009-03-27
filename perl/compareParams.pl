#!/usr/bin/perl -w

use DartSexpr;
use PhyloGram;

my $usage = "Usage: $0 <grammar files>\n\n";
die $usage unless @ARGV;

my @filename = @ARGV;
my (@gram, @param_hash);
for my $filename (@filename) {
    die $usage unless -e $filename;
    my $gram = PhyloGram->from_file ($filename);
    push @gram, $gram;
    my $param_hash = $gram->param_hash;
    push @param_hash, $param_hash;
#    warn $gram->grammar->name->value, " param_hash keys: ", join (" ", sort { $a cmp $b } keys %$param_hash);
}

my %param_hash_merge = map (%$_, @param_hash);
# sorting the parameters by the reverse of their name keeps things like "k" and "~k" together,
# as well as prefix-species constructs like "(. k_ NODE)"
my @params = sort { reverse($a) cmp reverse($b) } keys %param_hash_merge;

for my $param (@params) {
    my @info;
    push @info, "(param $param)";
    for (my $g = 0; $g < @gram; ++$g) {
	my $gram = $gram[$g];
	my @gramInfo;
	push @gramInfo, "(file " . $filename[$g] . ")";
	if (exists $param_hash[$g]->{$param}) {
	    push @gramInfo, "(def " . $param_hash[$g]->{$param}->value . ")";
	    my @obs = map ($_->find_all($param), $gram->grammar->find_all ("observed-counts"));
	    if (@obs) {
		die "More than one (observed-counts ...) expression found for $param" if @obs > 1;
		push @gramInfo, "(obs " . join ("/", $obs[0]->values) . ")";
	    }
	    push @info, "(@gramInfo)";
	}
    }
    print "(@info)\n";
}

