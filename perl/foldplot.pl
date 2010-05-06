#!/usr/bin/env perl -w

use GD;

# default opts
my ($xPre, $yPre, $xyPre, $xPost, $yPost, $xyPost, $xAxis, $yAxis);
my $outFile = "FOLDPLOT.png";
my $noText = 0;
my $dots = 0;
my $norm = 0;
my $useGray = 0;

# usage string
my $usage = "Usage: $0\n";
$usage   .= "        [-xpre <xPreFoldMatrix>]\n";
$usage   .= "        [-ypre <yPreFoldMatrix>]\n";
$usage   .= "        [-xypre <xyPreAlignMatrix>]\n";
$usage   .= "        [-xpost <xPostFoldMatrix>]\n";
$usage   .= "        [-ypost <yPostFoldMatrix>]\n";
$usage   .= "        [-xypost <xyPostAlignMatrix>]\n";
$usage   .= "        [-xtrue <xTrueFoldMatrix>]\n";
$usage   .= "        [-ytrue <yTrueFoldMatrix>]\n";
$usage   .= "        [-xytrue <xyTrueAlignMatrix>]\n";
$usage   .= "        [-out <outputFile>]    (default is $outFile)\n";
$usage   .= "        [-notext]   (doesn't label residues, creates a small image)\n";
$usage   .= "        [-dots]     (outlines \"true\" matches using dots, not solid lines)\n";
$usage   .= "        [-norm]     (normalises matrix before plotting)\n";
$usage   .= "        [-gray]     (uses grayscale, rather than heat)\n";
$usage   .= "\n";

# parse args
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-xpre") { $xPre = readMatrix (shift, \$xAxis, \$xAxis) }
    elsif ($arg eq "-ypre") { $yPre = readMatrix (shift, \$yAxis, \$yAxis) }
    elsif ($arg eq "-xypre") { $xyPre = readMatrix (shift, \$xAxis, \$yAxis) }
    elsif ($arg eq "-xpost") { $xPost = readMatrix (shift, \$xAxis, \$xAxis) }
    elsif ($arg eq "-ypost") { $yPost = readMatrix (shift, \$yAxis, \$yAxis) }
    elsif ($arg eq "-xypost") { $xyPost = readMatrix (shift, \$xAxis, \$yAxis) }
    elsif ($arg eq "-xtrue") { $xTrue = readMatrix (shift, \$xAxis, \$xAxis) }
    elsif ($arg eq "-ytrue") { $yTrue = readMatrix (shift, \$yAxis, \$yAxis) }
    elsif ($arg eq "-xytrue") { $xyTrue = readMatrix (shift, \$xAxis, \$yAxis) }
    elsif ($arg eq "-out") { $outFile = shift }
    elsif ($arg eq "-notext") { $noText = 1 }
    elsif ($arg eq "-dots") { $dots = 1 }
    elsif ($arg eq "-norm") { $norm = 1 }
    elsif ($arg eq "-gray") { $useGray = 1 }
    elsif ($arg =~ /^-/) { die $usage }
    else { push @argv, $arg }
}
@ARGV = @argv;
die $usage if @ARGV;
die $usage unless defined($xAxis) || defined($yAxis);

# get font
my $font = gdTinyFont;
my $charWidth = $font->width;
my $charHeight = $font->height;

# make it square
my ($cellWidth, $cellHeight) = ($charWidth, $charHeight);
if ($noText) { $cellWidth = $cellHeight = 2 }
else { if ($cellWidth > $cellHeight) { $cellHeight = $cellWidth } else { $cellWidth = $cellHeight } }

# get dimensions & some co-ords
my $lSize = (defined($yPre) || defined($yPost) || defined($xyPre)) ? (@$yAxis + 1) : 0;
my $rSize = (defined($xPre) || defined($xPost) || defined($xyPost)) ? (@$xAxis + 1) : 0;
my $uSize = (defined($yPre) || defined($yPost) || defined($xyPre)) ? (@$yAxis + 1) : 0;
my $dSize = (defined($xPre) || defined($xPost) || defined($xyPost)) ? (@$xAxis + 1) : 0;
my ($x0, $y0) = ($cellWidth, $cellHeight);
my ($x1, $y1) = ($cellWidth * $lSize, $cellHeight * $uSize);
my ($x2, $y2) = ($cellWidth * ($lSize + 1), $cellHeight * ($uSize + 1));
my ($x3, $y3) = ($cellWidth * ($lSize + $rSize), $cellHeight * ($uSize + $dSize));
my ($imWidth, $imHeight) = ($cellWidth * ($lSize + $rSize + 1), $cellHeight * ($uSize + $dSize + 1));

# create a new image
my $im = new GD::Image ($imWidth, $imHeight);
warn "Created image ($imWidth,$imHeight)\n";

# allocate some colors
my $black = $im->colorAllocate(0,0,0);
my $green = $im->colorAllocate(0,255,0);
my $blue = $im->colorAllocate(0,0,255);
my $yellow = $im->colorAllocate(180,180,0);
my $purple = $im->colorAllocate(180,0,180);
my $grey = $im->colorAllocate(16,16,16);
my $cyan = $im->colorAllocate(0,192,192);
my @col = ($black);
my $shades = 240;
for (my $i = 1; $i <= $shades; ++$i) {
    my $r = 32 + 224 * (3 * $i) / $shades;
    my $g = 256 * (3 * $i - $shades) / $shades;
    my $b = 256 * (3 * $i - 2*$shades) / $shades;
    if ($useGray) { $r = $g = $b = 255 * $i / ($shades + 1) }
    if ($r > 255) { $r = 255 } elsif ($r < 0) { $r = 0 }
    if ($g > 255) { $g = 255 } elsif ($g < 0) { $g = 0 }
    if ($b > 255) { $b = 255 } elsif ($b < 0) { $b = 0 }
    push @col, $im->colorAllocate ($r, $g, $b);
}
my $white = $col[@col-1];
my %charCol = ('a' => $purple, 't' => $blue, 'u' => $blue,
	       'c' => $green, 'g' => $yellow, '*' => $white);

# plot axes
if (defined($yPre) || defined($yPost)) { renderAxis ($im, $yAxis, $x0, $y0, $cellWidth, $cellHeight) }
if (defined($xPre) || defined($xPost)) { renderAxis ($im, $xAxis, $x2, $y2, $cellWidth, $cellHeight) }
if (defined($yPre)) { renderAxis ($im, $yAxis, 0, $y0, 0, $cellHeight) }
if (defined($yPost)) { renderAxis ($im, $yAxis, $x0, 0, $cellWidth, 0) }
if (defined($xPre)) { renderAxis ($im, $xAxis, $x2, $y3, $cellWidth, 0) }
if (defined($xPost)) { renderAxis ($im, $xAxis, $x3, $y2, 0, $cellHeight) }
if (defined($xyPre)) {
    renderAxis ($im, $xAxis, 0, $y2, 0, $cellHeight);
    renderAxis ($im, $yAxis, $x0, $y3, $cellWidth, 0);
}

if (defined($xyPost)) { renderAxis ($im, $xAxis, $x2, 0, $cellWidth, 0) }
if (defined($xyPost) || defined($xPost)) { renderAxis ($im, $xAxis, $x2, $y1, $cellWidth, 0) }
if (defined($xyPre) || defined ($xPre)) { renderAxis ($im, $xAxis, $x1, $y2, 0, $cellHeight) }
if (defined($yPost) || defined($xyPost)) { renderAxis ($im, $yAxis, $x1, $y0, 0, $cellHeight) }
if (defined($xyPost)) { renderAxis ($im, $yAxis, $x3, $y0, 0, $cellHeight) }
if (defined($xyPre) || defined($yPre)) { renderAxis ($im, $yAxis, $x0, $y1, $cellWidth, 0) }

# plot matrices
if (defined $xPre) { renderMatrix ($im, $xPre, $x2, $y2, 1, 1, $xTrue, $cyan) }
if (defined $yPre) { renderMatrix ($im, $yPre, $x0, $y0, 1, 1, $yTrue, $cyan) }
if (defined $xyPre) { renderMatrix ($im, $xyPre, $x0, $y2, 0, 1, $xyTrue, $cyan) }
if (defined $xPost) { renderMatrix ($im, $xPost, $x2, $y2, 1, 0, $xTrue, $cyan) }
if (defined $yPost) { renderMatrix ($im, $yPost, $x0, $y0, 1, 0, $yTrue, $cyan) }
if (defined $xyPost) { renderMatrix ($im, $xyPost, $x2, $y0, 0, 0, $xyTrue, $cyan) }

# Convert the image to PNG and print it to $outFile
local *OUTFILE;
open OUTFILE, ">$outFile" or die "Couldn't open '$outFile' for writing: $!";
binmode OUTFILE;
print OUTFILE $im->png;
close OUTFILE or die "Couldn't close '$outFile': $!";

# read matrix subroutine
sub readMatrix {
    my ($filename, $xAxisRef, $yAxisRef) = @_;

    local *FILE;
    open FILE, "<$filename" or die "Couldn't read post.prob. matrix from '$filename': $!";

    my $xAxis = <FILE>;
    chomp $xAxis;
    my ($dummy, @xAxis) = split /\s+/, $xAxis;

    if (defined $xAxisRef) {
	if (defined $$xAxisRef) {
	    my $oldAxis = $$xAxisRef;
	    die "Reading file '$filename': x axis length doesn't match" if @$oldAxis != @xAxis;
	    for (my $i = 0; $i < @xAxis; ++$i) {
		die "Reading file '$filename': x axis doesn't match" if $oldAxis->[$i] ne $xAxis[$i];
	    }
	} else {
	    $$xAxisRef = \@xAxis;
	}
    }

    my @data;
    my @yAxis;
    while (my $data = <FILE>) {
	chomp $data;
	my ($yCoord, @row) = split /\s+/, $data;
	push @yAxis, $yCoord;
	push @data, \@row;
	die "Reading file '$filename': row ", @yAxis+0, " has ", @row+0, " fields, expected ", @xAxis+0 if @row != @xAxis;
    }

    if (defined $yAxisRef) {
	if (defined $$yAxisRef) {
	    my $oldAxis = $$yAxisRef;
	    die "Reading file '$filename': y axis length doesn't match" if @$oldAxis != @yAxis;
	    for (my $i = 0; $i < @yAxis; ++$i) {
		die "Reading file '$filename': y axis doesn't match" if $oldAxis->[$i] ne $yAxis[$i];
	    }
	} else {
	    $$yAxisRef = \@yAxis;
	}
    }

    close FILE;

    if ($norm) {
	my ($min, $max);
	foreach my $row (@data) {
	    foreach my $datum (@$row) {
		$min = $datum if !defined($min) || $datum < $min;
		$max = $datum if !defined($max) || $datum > $max;
	    }
	}
	my $range = $max - $min;
	if ($range != 0) {
	    foreach my $row (@data) {
		foreach my $datum (@$row) {
		    $datum = ($datum - $min) / $range;
		}
	    }
	}
    }

    return \@data;
}

# subroutine to draw boxes, in outline, around nonzero cells of a matrix
# if diag != 0, render upper triangular only
# if flip != 0, then whole matrix is xy-inverted
sub renderOutlines {
    my ($im, $data, $xOrig, $yOrig, $upper, $flip, $col) = @_;
    if ($dots) { $im->setStyle ($col, $col, gdTransparent, gdTransparent) } else { $im->setStyle ($col) }
    for (my $y = 0; $y < @$data; ++$y) {
	my $row = $data->[$y];
	for (my $x = $upper ? $y+1 : 0; $x < @$row; ++$x) {
	    my $datum = $row->[$x];
	    if ($datum ne "." && $datum != 0) {
		my ($xpos, $ypos) = $flip ? ($y, $x) : ($x, $y);
		$im->rectangle ($xOrig + $xpos * $cellWidth, $yOrig + $ypos * $cellHeight,
				$xOrig + ($xpos + 1) * $cellWidth, $yOrig + ($ypos + 1) * $cellHeight,
				gdStyled);
	    }
	}
    }
}

# subroutine to render a matrix
# if upper != 0, render upper triangular only
# if flip != 0, then whole matrix is xy-inverted
# if true is defined, outline nonzero cells in color col
sub renderMatrix {
    my ($im, $data, $xOrig, $yOrig, $upper, $flip, $true, $col) = @_;
    warn "Drawing matrix from ($xOrig,$yOrig) to (", $xOrig+$cellWidth*@{$data->[0]}, ",", $yOrig+$cellHeight*@$data, ")\n";
    for (my $y = 0; $y < @$data; ++$y) {
	my $row = $data->[$y];
	for (my $x = $upper ? $y+1 : 0; $x < @$row; ++$x) {
	    my $datum = $row->[$x];
	    if ($datum eq ".") { $datum = 0 }
	    my $colIndex = $datum * @col;
	    if ($colIndex >= @col) { $colIndex = @col - 1}
	    elsif ($colIndex < 0) { $colIndex = 0 }
	    my ($xpos, $ypos) = $flip ? ($y, $x) : ($x, $y);
	    fillCell ($im, $xOrig + $xpos * $cellWidth, $yOrig + $ypos * $cellHeight, $col[$colIndex]);
	}
    }
    if (defined $true) {
	renderOutlines ($im, $true, $xOrig, $yOrig, $upper, $flip, $col);
    }
}

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
	    $im->char ($font, $x + $xOffset, $y + $yOffset, uc($c), $black);
	}
	$x += $xStep;
	$y += $yStep;
    }
}

# subroutine to render a diagonal in grey & white
sub renderDiag {
    my ($im, $x, $y, $cells) = @_;
    warn "Drawing diagonal from ($x,$y) to (", $x+$cellWidth*$cells-1, ",", $y+$cellHeight*$cells-1, ")\n";
    for (my $i = 0; $i < $cells; ++$i) {
	fillCell ($im, $x + $i * $cellWidth, $y + $i * $cellHeight, $grey);
    }
    $im->line ($x, $y, $x + $cells * $cellWidth - 1, $y + $cells * $cellHeight - 1, $white);
}

# low-level drawing subroutine to fill a cell
sub fillCell {
    my ($im, $x, $y, $col) = @_;
    $im->filledRectangle ($x, $y, $x + $cellWidth - 1, $y + $cellHeight - 1, $col);
}
