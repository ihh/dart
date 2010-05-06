#!/usr/bin/env perl -w

# set up pagers
my $less = "less -rF";
my $tail = "tail -f";
my $cat = "cat";

# get program name
my $progname = $0;
$progname =~ s/^.*\///;

# initialise help string & args
my $usage = "\n";
$usage .= "Usage:\n";
$usage .= "       $progname [options] <logfile(s)>\n";
$usage .= "   or\n";
$usage .= "       cat <logfile> | $progname [options]\n";
$usage .= "\n";
$usage .= "Options:\n";
$usage .= "\t          [-h]  this help message\n";
$usage .= "\t         [-fl]  print file & line info\n";
$usage .= "\t         [-dt]  print date & time info\n";
$usage .= "\t        [-lev]  print log level\n";
$usage .= "\t        [-tag]  print tags\n";
$usage .= "\t   [-f,-fancy]  synonym for \"-fl -dt -lev -tag\"\n";
$usage .= "\t        [-cat]  pipe file through '$cat' instead of '$less'\n";
$usage .= "\t       [-tail]  pipe file from '$tail' and through '$cat' instead of '$less'\n";
$usage .= "\t[-min <level>]  only display log messages above this log level\n";
$usage .= "\t [-skip <tag>]  skip log messages with this tag\n";
$usage .= "\t [-show <tag>]  show log messages with this tag (overrides -min and -skip)\n";

my $print_file_and_line = 0;
my $print_date_and_time = 0;
my $print_level = 0;
my $print_tags = 0;
my $pipe_from_tail = 0;
my $pipe_to_cat = 0;
my $tab = "";
my $minlev;
my %show;
my %skip;

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-lev") {
	$print_level = 1;
    } elsif ($arg eq "-fl") {
	$print_file_and_line = 1;
    } elsif ($arg eq "-dt") {
	$print_date_and_time = 1;
    } elsif ($arg eq "-tag") {
	$print_tags = 1;
    } elsif ($arg eq "-f" || $arg eq "-fancy") {
	unshift @ARGV, qw(-fl -dt -lev -tag);
    } elsif ($arg eq "-cat") {
	$pipe_to_cat = 1;
    } elsif ($arg eq "-tail") {
	$pipe_from_tail = $pipe_to_cat = 1;
    } elsif ($arg eq "-min") {
	defined ($minlev = shift) or die $usage;
    } elsif ($arg eq "-skip") {
	my $tag = shift;
	defined ($tag) or die $usage;
	$skip{$tag} = 1;
    } elsif ($arg eq "-show") {
	my $tag = shift;
	defined ($tag) or die $usage;
	$show{$tag} = 1;
    } else {
	push @argv, $arg;
    }
}
@argv = qw(-) unless @argv;

# open less pager
my $pager = $pipe_to_cat ? $cat : $less;
local *LESS;
open (LESS, "|$pager") or open (LESS, ">-");

# autoflush
select LESS;
$| = 1;
select STDOUT;

# set up ANSI color codes
my @col = (map (e(7).e(23).e($_), reverse (31..35), 37),  # inverse, not-italics, (magenta..red,white)
	   map (e(27).e(23).e($_), 31..37),  # not-inverse, not-italics, (red..white)
	   e(7).e(23).e(36));  # inverse, not-italics, cyan


my $w = e(27).e(23).e(37);  # not-inverse, not-italics, white
my $invw = e(7).e(23).e(37);  # inverse, not-italics, white
my $invr = e(7).e(23).e(31);  # inverse, not-italics, red
my $invg = e(7).e(23).e(32);  # inverse, not-italics, green
my $invy = e(7).e(23).e(33);  # inverse, not-italics, yellow
my $invb = e(7).e(23).e(34);  # inverse, not-italics, blue
my $invm = e(7).e(23).e(35);  # inverse, not-italics, magenta
my $invc = e(7).e(23).e(36);  # inverse, not-italics, cyan

# parse logfile
my ($date, $time, $file, $line, $tags) = map ("", 1..5);
my ($last_date, $last_time, $last_file, $last_line, $last_tags) = map ("", 1..5);
my ($lev, $last_lev) = (9, 9);
my @tag;
my $c = "";
foreach my $infile (@argv) {
    local *INFILE;
    if ($pipe_from_tail) {
	open INFILE, "$tail $infile|" or die "Couldn't run '$tail $infile'";
    } else {
	open INFILE, "<$infile" or die "Couldn't open $infile";
    }
    while (<INFILE>) {
	my $text = $_;
	if ($text =~ /<log.*>/) {
	    if (/date="(\S+)"/) { $date = $1 }
	    if (/time="(\S+)"/) { $time = $1 }
	    if (/file="(\S+)"/) { $file = $1 }
	    if (/line=(\d+)/)   { $line = $1 }
	    if (/level=(\-?\d+)/) { $lev = $1 }

	    $c = $col [$lev<-5 ? 0 : ($lev+5>=@col ? @col-1 : $lev+5)];
	    if (/tags=\"([^\"]+)\"/) {
		$tags = $1;
		@tag = split /\s+/, $tags;
	    } else {
		$tags = "";
		@tag = ();
	    }
	} else {
	    my $show = grep (exists($show{$_}), @tag);
	    next if !$show && defined($minlev) && $lev < $minlev;
	    next if !$show && %skip && grep (exists($skip{$_}), @tag);
	    my $info = 0;
	    if ($print_level) {	print LESS $invc, ($lev == $last_lev ? ("  ") : (($lev<0||$lev>9?():(' ')), $lev)), $invc, $w, ' ' }
	    if ($print_file_and_line && ($file ne $last_file || $line ne $last_line)) { print LESS $invr, $file, '#', $line, $invr, $w; $last_file = $file; $last_line = $line; $info = 1 }
	    if ($print_date_and_time && $time ne $last_time) { print LESS $invg, $time, $invg, $w; $last_time = $time; $info = 1 }
	    if ($print_date_and_time && $date ne $last_date) { print LESS $invm, $date, $invm, $w; $last_date = $date; $info = 1 }
	    if ($print_tags && $tags ne $last_tags) { print LESS $invy, $tags, $invy, $w; $last_tags = $tags; $info = 1 }
	    print LESS "\t" if $info;
	    chomp $text;
	    print LESS $w, $c, $text, $c, $w, "\n";  # sandwiching with $w$c allows paging back with 'less'
	}
    }
    close INFILE;
}

# escape code subroutine
sub e {
    my $code = shift;
    return chr(27)."[$code"."m";
}
