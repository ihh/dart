#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use GD;

# usage string
my $usage = "Usage: $0 [-diag] [-noscale] [-notext] <file>\n";

# default opts
my $diag = 0;
my $printScale = 1;
my $noText = 0;

# parse args
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-diag") { $diag = 1 }
    elsif ($arg eq "-noscale") { $printScale = 0 }
    elsif ($arg eq "-notext") { $noText = 1 }
    elsif ($arg =~ /^-/) { die $usage }
    else { push @argv, $arg }
}
@ARGV = @argv;

# read data
my $xAxis = <>;
chomp $xAxis;
my ($dummy, @xAxis) = split /\s+/, $xAxis;

my @data;
my @yAxis;
while (my $data = <>) {
    chomp $data;
    my ($yCoord, @row) = split /\s+/, $data;
    push @yAxis, $yCoord;
    push @data, \@row;
}

# if -diag specified, check dimensions match
if ($diag && @xAxis != @yAxis) { die "-diag option specified, but x & y sequences don't match" }

# get font
my $font = gdTinyFont;
my $charWidth = $font->width;
my $charHeight = $font->height;

# make it square
my ($cellWidth, $cellHeight) = ($charWidth, $charHeight);
if ($noText) { $cellWidth = $cellHeight = 2 }
else { if ($cellWidth > $cellHeight) { $cellHeight = $cellWidth } else { $cellWidth = $cellHeight } }

# create a new image
my $imWidth = $cellWidth * (@xAxis + 2);
my $imHeight = $cellHeight * (@yAxis + ($printScale ? 3 : 2));
my $im = new GD::Image ($imWidth, $imHeight);

# allocate some colors
my $black = $im->colorAllocate(0,0,0);
my $green = $im->colorAllocate(0,255,0);
my $blue = $im->colorAllocate(0,0,255);
my $yellow = $im->colorAllocate(180,180,0);
my $purple = $im->colorAllocate(180,0,180);
my $grey = $im->colorAllocate(16,16,16);
my @col = ($black);
my $shades = 250;
for (my $i = 1; $i <= $shades; ++$i) {
    my $r = 32 + 224 * (3 * $i) / $shades;
    my $g = 256 * (3 * $i - $shades) / $shades;
    my $b = 256 * (3 * $i - 2*$shades) / $shades;
    if ($r > 255) { $r = 255 } elsif ($r < 0) { $r = 0 }
    if ($g > 255) { $g = 255 } elsif ($g < 0) { $g = 0 }
    if ($b > 255) { $b = 255 } elsif ($b < 0) { $b = 0 }
    push @col, $im->colorAllocate ($r, $g, $b);
}
my $white = $col[@col-1];
my %charCol = ('a' => $purple, 't' => $blue, 'u' => $blue,
	       'c' => $green, 'g' => $yellow, '*' => $white);

# plot color scale
if ($printScale) {
    my $ltext = "Probability: 0 ";
    my $rtext = " 1";
    my $lchars = length ($ltext);
    my $rchars = length ($rtext);
    my $barWidth = ($imWidth - $charWidth * ($lchars + $rchars)) / @col;
    for (my $c = 0; $c < @col; ++$c) {
	$im->filledRectangle ($charWidth * $lchars + $c * $barWidth, $imHeight - $charHeight,
			      $charWidth * $lchars + ($c + 1) * $barWidth - 1, $imHeight,
			      $col[$c]);
    }
    $im->string ($font, 0, $imHeight - $charHeight, $ltext, $white);
    $im->string ($font, $imWidth - $charWidth * $rchars, $imHeight - $charHeight, $rtext, $white);
}


# plot data
for (my $y = 0; $y < @yAxis; ++$y) {
    for (my $x = $diag ? $y+1 : 0; $x < @xAxis; ++$x) {
	my $datum = $data[$y]->[$x];
	if ($datum eq ".") { $datum = 0 }
	my $colIndex = int ($datum * @col);
	if ($colIndex >= @col) { $colIndex = @col - 1}
	elsif ($colIndex < 0) { $colIndex = 0 }
	$im->filledRectangle (($x + 1) * $cellWidth, ($y + 1) * $cellHeight,
			      ($x + 2) * $cellWidth - 1, ($y + 2) * $cellHeight - 1,
			      $col[$colIndex]);
    }
}

# plot sequence axes
if ($diag) {
    my $len = @xAxis;
#    my $triangle = new GD::Polygon;
#    $triangle->addPt (0, 0);
#    $triangle->addPt (0, ($len+1)*$cellHeight);
#    $triangle->addPt (($len+1)*$cellWidth, ($len+1)*$cellHeight);
#    $im->filledPolygon ($triangle, $white);

    renderAxis ($im, \@xAxis, $cellWidth, 0, $cellWidth, 0);
    renderAxis ($im, \@xAxis, $cellWidth, $cellHeight, $cellWidth, $cellHeight);
    renderAxis ($im, \@xAxis, $imWidth - $cellWidth, $cellHeight, 0, $cellHeight);
} else {
    renderAxis ($im, \@xAxis, $cellWidth, 0, $cellWidth, 0);
    renderAxis ($im, \@xAxis, $cellWidth, $imHeight - ($printScale ? 2 : 1) * $cellHeight, $cellWidth, 0);
    renderAxis ($im, \@yAxis, 0, $cellHeight, 0, $cellHeight);
    renderAxis ($im, \@yAxis, $imWidth - $cellWidth, $cellHeight, 0, $cellHeight);
}

# make sure we are writing to a binary stream
binmode STDOUT;

# Convert the image to PNG and print it on standard output
print $im->png;

# subroutine to render an axis
sub renderAxis {
    my ($im, $axisRef, $x, $y, $xStep, $yStep) = @_;
    warn "Drawing sequence from ($x,$y) to (", $x+$xStep*@$axisRef, ",", $y+$yStep*@$axisRef, ")\n";
    my $xOffset = ($cellWidth - $charWidth) / 2;
    my $yOffset = ($cellHeight - $charHeight) / 2;
    for (my $i = 0; $i < @$axisRef; ++$i) {
	my $c = $axisRef->[$i];
	fillCell ($im, $x, $y, $charCol{lc$c});
	if ($cellWidth >= $charWidth && $cellHeight >= $charHeight) {
	    $im->char ($font, $x + $xOffset, $y + $yOffset, $c, $black);
	}
	$x += $xStep;
	$y += $yStep;
    }
}

# low-level drawing subroutine to fill a cell
sub fillCell {
    my ($im, $x, $y, $col) = @_;
    $im->filledRectangle ($x, $y, $x + $cellWidth - 1, $y + $cellHeight - 1, $col);
}
