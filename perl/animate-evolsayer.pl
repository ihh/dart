#!/usr/bin/perl -w

# options
my $frames_per_sub = 50;
my $mpegFile = "evolsayer.mpg";
my $paramFile = "PARAM";  # parameter file that will be generated for mpeg_encode
my $frameDir = "frames";  # frame directory
my $convert_size = "200x200";  # size of frame images
my $savePPM = 0;
my $rnaplot_type = 1;
my $tkfst = 1;

# file suffices
my $imgSuffix = "_ss.ps";  # suffix added by RNAplot
my $ppmSuffix = "_ss.ppm";  # suffix of converted PPM file

# helpers
my $rnaplot = "RNAplot";
my $convert = "convert";
my $mpeg_encode = "mpeg_encode";

# usage message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "$progname: animation of evolsayer simulation .history files\n";
$usage .=   "\n";
$usage .=   "Usage: $progname <.history file>\n";
$usage .=   "              [-h]  print this message\n";
$usage .=   "           [-rate]  frames per unit of evolutionary time (default is $frames_per_sub)\n";
$usage .=   "        [-notkfst]  don't show TKFST string representation\n";
$usage .=   "    [-mpeg <file>]  MPEG filename (default is $mpegFile)\n";
$usage .=   " [-rnaplot <path>]  path to Vienna RNAplot program (default is '$rnaplot')\n";
$usage .=   " [-convert <path>]  path to ImageMagick PNG-->PPM conversion program (default is '$convert')\n";
$usage .=   "  [-encode <path>]  path to Berkeley MPEG Encoder (default is '$mpeg_encode')\n";
$usage .=   "         [-radial]  use radial layout in RNAplot, instead of naview layout\n";
$usage .=   "   [-size <w>x<h>]  size of movie frame images (default is $convert_size)\n";
$usage .=   "   [-param <file>]  generated parameter file for Berkeley MPEG Encoder (default is '$paramFile')\n";
$usage .=   " [-framedir <dir>]  directory in which to store frame images (default is '$frameDir')\n";
$usage .=   "        [-saveppm]  don't delete PPM frame image files\n";
$usage .=   "\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-rate") {
	defined ($frames_per_sub = shift) or die $usage;
    } elsif ($arg eq "-notkfst") {
	$tkfst = 0;
    } elsif ($arg eq '-mpeg') {
	defined ($mpegFile = shift) or die $usage;
    } elsif ($arg eq '-rnaplot') {
	defined ($rnaplot = shift) or die $usage;
    } elsif ($arg eq '-convert') {
	defined ($convert = shift) or die $usage;
    } elsif ($arg eq '-encode') {
	defined ($mpeg_encode = shift) or die $usage;
    } elsif ($arg eq '-radial') {
	$rnaplot_type = 0;
    } elsif ($arg eq '-size') {
	defined ($convert_size = shift) or die $usage;
    } elsif ($arg eq '-param') {
	defined ($paramFile = shift) or die $usage;
    } elsif ($arg eq '-framedir') {
	defined ($frameDir = shift) or die $usage;
    } elsif ($arg eq '-saveppm') {
	$savePPM = 1;
    } else {
	push @argv, $arg;
    }
}
unless (@argv) {
    @argv = ('-');
    warn "[waiting for .history file on standard input]\n";
}
die $usage unless @argv == 1;
my ($filename) = @argv;

# let's make a movie
mkdir ($frameDir) unless -d $frameDir;

# open history file
local *HISTORY;
open HISTORY, "<$filename" or die $!;

# write frames
my @ppm;
my $last_time = 0;
my $n_frame = 0;
while (<HISTORY>) {
    my ($time, $tkfst, $seq, $ss) = split;
    if (@ppm && $last_time + 1/$frames_per_sub < $time) {
	my $n_frames = ($time - $last_time) * $frames_per_sub;
	if ($n_frames >= 1) {
	    push @ppm, map ($ppm[-1], 1..$n_frames);
	}
    }
    $last_time = $time;

    $ss =~ tr/<>/()/;  # convert Stockholm SS_cons line to the format RNAplot likes
    
    my $prefix = "$frameDir/$n_frame";
    warn "Drawing frame at time $time ($prefix)\n";

    my $rnaplot_command = "$rnaplot -o ps -t $rnaplot_type";
    $rnaplot_command .= " --post '" . centered_text_ps($tkfst) . "'" if $tkfst;

    local *RNAPLOT;
    open RNAPLOT, "| $rnaplot_command";

    for my $line (">$prefix", $seq, $ss, '@') {
	print RNAPLOT "$line\n";
#	warn "[rnaplot] $line\n";
    }

    close RNAPLOT;

    my $img = "$prefix$imgSuffix";
    my $ppm = "$prefix$ppmSuffix";

    system "$convert -size $convert_size $img $ppm";
    push @ppm, $ppm;

    ++$n_frame;
}

# run mpeg_encode
mpegEncode ($mpegFile, \@ppm);

# delete PPM's
warn "(deleting PPM image files)\n";
unless ($savePPM) {
    for my $ppm (@ppm) {
	unlink $ppm;
    }
}

# mpegEncode
sub mpegEncode {
    my ($mpegFile, $frameListRef) = @_;
    local *PARAM;
    open PARAM, ">$paramFile" or die "Couldn't open param file '$paramFile': $!";

    print PARAM "PATTERN          IBBPBBPBBPBBPBBP\n";
    print PARAM "OUTPUT           $mpegFile\n";
    print PARAM "BASE_FILE_FORMAT PNM\n";
    print PARAM "INPUT_CONVERT    *\n";
    print PARAM "GOP_SIZE         16\n";
    print PARAM "SLICES_PER_FRAME 1\n";

    print PARAM "INPUT_DIR        .\n";
    print PARAM "INPUT\n";
    print PARAM map ("$_\n", @$frameListRef);
    print PARAM "END_INPUT\n";

    print PARAM "PIXEL            HALF\n";
    print PARAM "RANGE            10\n";
    print PARAM "PSEARCH_ALG      LOGARITHMIC\n";
    print PARAM "BSEARCH_ALG      CROSS2\n";
    print PARAM "IQSCALE          8\n";
    print PARAM "PQSCALE          10\n";
    print PARAM "BQSCALE          25\n";
    print PARAM "REFERENCE_FRAME  ORIGINAL\n";
    print PARAM "BIT_RATE         1000000\n";   # constant bit-rate
#    print PARAM "BUFFER_SIZE      327680\n";  # buffer size
    print PARAM "FRAME_RATE       30\n";   # frames per second
    print PARAM "FORCE_ENCODE_LAST_FRAME\n";

    close PARAM or die "Couldn't close param file '$paramFile': $!";

    system "$mpeg_encode $paramFile";
}

# Postscript for drawing TKFST string
# The following are defined in the RNAplot-generated Postscript:
#  RED: { 1 0 0 }
#  fsize: font size
#  cshow: show string on top of stack, centered (uses fsize)
sub centered_text_ps {
    my ($text) = @_;
    return join (" ",
		 'RED setrgbcolor',
		 '/fsize size', "($text)", 'length div 1.5 mul def',
		 '/Helvetica findfont fsize scalefont setfont',
		 'xmin xmax add 2 div',
		 'ymin ymax add 2 div',
		 'moveto', "($text)", 'cshow');
}
