#!/usr/bin/perl -w

my $bindir = "../../../bin";
my $xrate = "$bindir/xrate";

my $failSuffix = ".fail";

my $dummyAlign = "dummy.stock";

local *DIR;
opendir DIR, '.';
my @failGramFile = grep -f $_ && /$failSuffix$/, readdir DIR;
closedir DIR;

my ($npass, $ntot) = (0, 0);
for my $failGramFile (@failGramFile) {
    syswarn ("$xrate -g $failGramFile $dummyAlign -noa >/dev/null");
    if ($? == 0) {
	warn "\nExpected grammar '$failGramFile' to fail, but it didn't...\n\n";
    } else {
	++$npass;
    }
    ++$ntot;
}
warn "$npass/$ntot grammars failed on cue\n";
print $npass==$ntot ? "ok\n" : "not ok\n";

sub syswarn {
    my ($command) = @_;
    warn "\n***** Running '$command' *****\n";
    return `$command`;
}
