#!/usr/bin/perl

# Strip off the "noself-" prefix from the filename of a Perl grammar generating script,
# add the command-line option "-noself", and run on the default template file for that script.

my $tmplSuffix = ".tmpl.eg";
$0 =~ /^(.*)noself-([^\/]+)(\.pl)$/
    or die "Program name must be of the form PATH/noself-PROGNAME.pl";
my ($path, $progName, $plSuffix) = ($1, $2, $3);
my $command = "perl $path$progName$plSuffix -noself $path$progName$tmplSuffix";
warn $command, "\n";
exec $command;

