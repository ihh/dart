#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm::Database;
use Newick;
use POSIX qw(ceil floor);

# options
my $delay_secs = .02;  # delay in seconds between "frames"
my $cols;  # alignment columns per row
my $pixels;  # use pixels, not text
my $mpegFile;
my $paramFile = "PARAM";  # parameter file that will be generated for mpeg_encode
my $frameDir = "frames";  # frame directory
my $imagePrefix = "img";  # frame image filename prefix
my $savePPM = 0;
my $useSVG = 0;

# Rasmol color scheme
my %rasRGB = ('D' => [230,10,10],
	      'E' => [230,10,10],
	      'C' => [230,230,0],
	      'M' => [230,230,0],
	      'K' => [20,90,255],
	      'R' => [20,90,255],
	      'S' => [250,150,0],
	      'T' => [250,150,0],
	      'F' => [50,50,170],
	      'Y' => [50,50,170],
	      'N' => [0,220,220],
	      'Q' => [0,220,220],
	      'G' => [235,235,235],
	      'L' => [15,130,15],
	      'V' => [15,130,15],
	      'I' => [15,130,15],
	      'A' => [200,200,200],
	      'W' => [180,90,180],
	      'H' => [130,130,210],
	      'P' => [220,150,130],
	      'B' => [160,160,160],
	      'Z' => [160,160,160],
	      'X' => [140,140,140],
	      '*' => [60,60,60]);

# file suffices
my $pngSuffix = ".png";
my $ppmSuffix = ".ppm";
my $svgSuffix = ".svg";

# helpers
my $clear = "clear";
my $convert = "convert";
my $mpeg_encode = "mpeg_encode";

# help message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "$progname: animation of Stockholm alignment traces from MCMC sampling runs\n";
$usage .=   "\n";
$usage .=   "Usage: $progname <Stockholm file>\n";
$usage .=   "              [-h]  print this message\n";
$usage .=   "          [-delay]  frame delay in seconds (default is $delay_secs)\n";
$usage .=   "       [-cols <N>]  number of columns per row\n";
$usage .=   "         [-pixels]  draw characters as pixels, not text\n";
$usage .=   "    [-mpeg <file>]  write MPEG to file\n";
$usage .=   "   [-clear <path>]  path to 'clear screen' program (default is '$clear')\n";
$usage .=   " [-convert <path>]  path to ImageMagick PNG-->PPM conversion program (default is '$convert')\n";
$usage .=   "  [-encode <path>]  path to Berkeley MPEG Encoder (default is '$mpeg_encode')\n";
$usage .=   "   [-param <file>]  generated parameter file for Berkeley MPEG Encoder (default is '$paramFile')\n";
$usage .=   " [-framedir <dir>]  directory in which to store frame images (default is '$frameDir')\n";
$usage .=   "        [-saveppm]  don't delete PPM frame image files\n";
$usage .=   "            [-svg]  use SVG instead of PNG format for frame image files\n";
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
    } elsif ($arg eq '-pixels') {
	$pixels = 1;
    } elsif ($arg eq '-mpeg') {
	defined ($mpegFile = shift) or die $usage;
    } elsif ($arg eq '-clear') {
	defined ($clear = shift) or die $usage;
    } elsif ($arg eq '-convert') {
	defined ($convert = shift) or die $usage;
    } elsif ($arg eq '-encode') {
	defined ($mpeg_encode = shift) or die $usage;
    } elsif ($arg eq '-param') {
	defined ($paramFile = shift) or die $usage;
    } elsif ($arg eq '-framedir') {
	defined ($frameDir = shift) or die $usage;
    } elsif ($arg eq '-saveppm') {
	$savePPM = 1;
    } elsif ($arg eq '-svg') {
	$useSVG = 1;
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

# do we need GD or GD::SVG?
my ($imagePackage, $imageClass, $fontClass, $fontMethod, $imageOutputMethod, $imageFileSuffix);
if (defined $mpegFile) {
    $imagePackage = $useSVG ? "GD::SVG" : "GD";
    $imageClass = $imagePackage . "::Image";

    # Changed $fontClass from "GD::Font" to "GD" since this now seems (empirically) to be the norm... in stark contrast to the documentation! IH, 4/7/09
    # $fontClass = $imagePackage . "::Font";
    # $fontMethod = "Giant";
    $fontClass = $imagePackage;
    $fontMethod = "gdGiantFont";

    $imageOutputMethod = $useSVG ? "svg" : "png";
    $imageFileSuffix = $useSVG ? $svgSuffix : $pngSuffix;

    eval "use $imagePackage";
}

# read alignments
my $db = Stockholm::Database->from_file ($filename);

# let's make a movie
my $ticks_per_sec;  # make sure this is globally defined
if (defined $mpegFile) {
    # mpegs
    make_mpeg();
} else {
    # text animations
    die "If -pixels is specified, you must also specify the -mpeg option" if $pixels;
    warn "Animation has ", @$db+0, " alignment frames & will take ", $delay_secs * @$db, " seconds (plus display overhead)\n";
    warn "Calibrating delay loop... (should take 1 second)\n";
    $ticks_per_sec = 0;
    alarm 1;
    $SIG{'ALRM'} = \&alarm_sub;   # hmm, must be a less hacky way of doing this...
    while (1) { ++$ticks_per_sec }   # alarm should interrupt this infinite loop and send us to alarm_sub()
}

# called by text animation (via 'alarm')
sub alarm_sub {
    alarm 0;  # turn off that alarm
    my $delay_ticks = $delay_secs * $ticks_per_sec;

    for my $stock (@$db) {
	my $next_frame = $stock->to_string ($cols);
	system $clear;  # clear the screen
	print $next_frame;
	for (my $ticks = 0; $ticks < $delay_ticks; ++$ticks) { }
    }

    exit;  # this is the hackiest part. without this 'exit', seems to hang forever (in the infinite delay loop?)
}

# MPEGs
sub make_mpeg {
    mkdir ($frameDir) unless -d $frameDir;

    # prepare text for all frames
    my @frame;
    my ($xsize, $ysize);
    my ($cw, $ch, $font);
    my $columns = max (map ($_->columns, @$db));  # columns
    my $rows = max (map (@{$_->seqname} + 0, @$db));  # rows
    if ($pixels) {
	# compact: 1 pixel per alignment char
	$cols = ceil (sqrt ($columns * $rows)) unless defined $cols;  # try to make image square
	$xsize = min ($cols, $columns);
	$ysize = $rows * int ($columns / $cols + 1);
	$columns = ceil ($columns / $cols) * $cols;  # round up $columns to use all the space
    } else {
	# expanded: ($font->width * $font->height) pixels per alignment char
	my ($width, $height) = (0, 0);
	for my $stock (@$db) {
	    warn "Preparing frame ", @frame+0, "\n";
	    my $stock_text = $stock->to_string ($cols);
	    my @stock_lines = split /\n/, $stock_text;
	    push @frame, \@stock_lines;
	    $width = max ($width, map (length, @stock_lines));
	    $height = max ($height, @stock_lines + 0);
	}

	# font & dimensions
	$font = $fontClass->$fontMethod;
	$cw = $font->width;
	$ch = $font->height;
	($xsize, $ysize) = ($width * $cw, $height * $ch);
    }

    # image & colors
    my $im = $imageClass->new ($xsize, $ysize);
    my $black = $im->colorAllocate (0, 0, 0);
    my $white = $im->colorAllocate (255, 255, 255);

    # Rasmol colors
    my %imColor = map (($_ => $im->colorAllocate (@{$rasRGB{$_}})), keys %rasRGB);

    # write frames
    my @ppm;
    for (my $n = 0; $n < @$db; ++$n) {
	warn "Drawing frame $n/", @$db-1, "\n";

	clear ($im, $xsize, $ysize, $black);

	if ($pixels) {
	    my $stock = $db->[$n];
	    my @seqname = sort @{$stock->seqname};
	    my @rowOrder = (0 .. @seqname - 1);
	    if (@{$stock->gf_NH}) {
		my $tree = Newick->from_string (join ("", @{$stock->gf_NH}));
		my %rowIndex = map (($seqname[$_] => $_), @rowOrder);
		@rowOrder = map ((defined($_) ? ($_) : ()), map ($rowIndex{$_}, @{$tree->node_name}));
	    }
	    for my $row (@rowOrder) {
		my $s = $stock->seqdata->{$seqname[$row]};
		my $mul = $columns / $stock->columns;
		for my $col (0 .. $stock->columns - 1) {
		    my $c = uc substr ($s, $col, 1);
		    if (exists $imColor{$c}) {
			for my $x (int ($mul * $col) .. int ($mul * ($col + 1)) - 1) {
			    $im->setPixel ($x % $cols, $row + $rows * int ($x / $cols), $imColor{$c});
			}
		    }
		}
	    }
	} else {
	    for (my $line = 0; $line < @{$frame[$n]}; ++$line) {
		drawString ($im, $font, 0, $line * $ch, $cw, $frame[$n]->[$line], $white);
	    }
	}

	my $prefix = "$frameDir/$imagePrefix$n";

	my $img = "$prefix$imageFileSuffix";
	my $ppm = "$prefix$ppmSuffix";

	local *IMG;
	open IMG, ">$img" or die "Couldn't write to '$img': $!";
	binmode IMG;
	print IMG $im->$imageOutputMethod;
	close IMG or die "Couldn't close '$img': $!";

	system "$convert -depth 8 $img $ppm";
	push @ppm, $ppm;
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
}

# min
sub min {
    my ($min, @x) = @_;
    for my $x (@x) { $min = $x if $x < $min }
    return $min;
}

# max
sub max {
    my ($max, @x) = @_;
    for my $x (@x) { $max = $x if $x > $max }
    return $max;
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

# subroutine to draw a string
# in PNG mode, just calls GD::Image::string once
# in SVG mode, calls it one character at a time (somehow SVG can't seem to manage fixed-width fonts...)
sub drawString {
    my ($im, $font, $x, $y, $cw, $string, $color) = @_;
    if ($useSVG) {
	for (my $i = 0; $i < length($string); ++$i) {
	    my $c = substr ($string, $i, 1);
	    if ($c ne " ") {
		$im->string ($font, $x + $i * $cw, $y, $c, $color);
	    }
	}
    } else {
	$im->string ($font, $x, $y, $string, $color);
    }
}

# subroutine to blank the frame
sub clear {
    my ($im, $xsize, $ysize, $bgcol) = @_;
    if ($useSVG) {
	allocateImage();
    } else {
	$im->filledRectangle (0, 0, $xsize, $ysize, $bgcol);
    }
}


