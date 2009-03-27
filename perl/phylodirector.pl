#!/usr/bin/perl -w

# prepend script dir to path
my $scriptDir = $0;
$scriptDir =~ s!/[^/]+$!!;
push @INC, $scriptDir;

use Stockholm;
use Newick;

# constants
my $pi = 2 * atan2(1,0);  # Pi = 3.14159...

# helpers
my $convert = "convert";
my $mpeg_encode = "mpeg_encode";

# file suffices
my $pngSuffix = ".png";
my $ppmSuffix = ".ppm";
my $svgSuffix = ".svg";

# directories, filenames
my $colImagePrefix = "col";
my $paramFile = "PARAM";
my $frameDir = "frames";
my $mpegFile = "ehmm.mpg";
my $collapsedMpegFile = "collapsed-ehmm.mpg";

# other params
my $savePPM = 0;
my $useSVG = 0;
my $verbose = 0;

# overall scale
my $scale = 1;

# speed factors
my $flash = 1;
my $pauseTime = 5;
my $emitPauseFlashRatio = 1;
my $windPauseFlashRatio = 2;

# EHMM icon dimensions
my $bxsize = 100 * $scale;  # branch x-size
my $bysize = 50 * $scale;  # branch y-size
my $ssize = 20 * $scale;  # state half-size
my $nsize = 5 * $scale;  # node radius
my $bwidth = 1 * $scale;  # branch half-width
my $alen = 20 * $scale;  # branch arrow length
my $awidth = 8 * $scale;  # branch arrow half-width
my $anaoff = 5 * $scale;  # active-node arrow offset
my $anllen = 20 * $scale;  # active-node line length
my $anlwidth = 2 * $scale;  # active-node line half-width
my $analen = 15 * $scale;  # active-node arrow length
my $anawidth = 10 * $scale;  # active-node arrow half-width

# pixel size multiplier for mpeg_encode
my $mpegPixelMul = 16;

# flag to leave space for EHMM zero-node
my $useZeroSpace = 0;

# start, end columns
my $startCol = 1;
my $endCol;

# flags specifying whether to build animations
my $buildExpanded = 1;
my $buildMpegs = 1;

# usage text
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage .= "Usage: $progname [options] [Stockholm alignment file]\n";
$usage .= "\n";
$usage .= "   -convert <path>  path to ImageMagick PNG-->PPM conversion program (default is '$convert')\n";
$usage .= "    -encode <path>  path to Berkeley MPEG Encoder (default is '$mpeg_encode')\n";
$usage .= "   -framedir <dir>  directory in which to store frame images (default is '$frameDir')\n";
$usage .= "     -param <file>  parameter file for Berkeley MPEG Encoder (default is '$paramFile')\n";
$usage .= "  -expanded <file>  MPEG filename for 'expanded' animation (default is '$mpegFile')\n";
$usage .= " -collapsed <file>  MPEG filename for 'collapsed' animation (default is '$collapsedMpegFile')\n";
$usage .= "         -noexpand  don't generate 'expanded' MPEG\n";
$usage .= "           -nompeg  don't generate MPEGs, just make still frames\n";
$usage .= "      -start <col>  start column (default is $startCol)\n";
$usage .= "        -end <col>  end column (default is last column of alignment)\n";
$usage .= "       -fflash <n>  number of frames in each 'on' and 'off' state of a flash (default is $flash)\n";
$usage .= "       -nflash <n>  base number of flashes (default is $pauseTime)\n";
$usage .= "         -emit <n>  multiplier for number of flashes during emission cascades (default is $emitPauseFlashRatio)\n";
$usage .= "         -wind <n>  multiplier for number of flashes during transitions (default is $windPauseFlashRatio)\n";
$usage .= "        -scale <n>  size scaling factor\n";
$usage .= "          -saveppm  don't delete PPM frame image files\n";
$usage .= "              -svg  use SVG instead of PNG format for frame image files\n";
$usage .= "          -verbose  print lots of stuff to standard error\n";
$usage .= "\n";
$usage .= "Notes:\n";
$usage .= " [1] Stockholm alignment must contain a Newick-format tree, prefixed by tag '#=GF NH'\n";
$usage .= " [2] All nodes of the tree will be named (including internal nodes)\n";
$usage .= " [3] Every node name should correspond to the name of a sequence in the alignment\n";
$usage .= " [4] All internal tree nodes must be binary except for the root node, which can be unary\n";
$usage .= "\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h" || $arg eq "-help" || $arg eq "--help") {
	die $usage;
    } elsif ($arg eq '-convert') {
	defined ($convert = shift) or die $usage;
    } elsif ($arg eq '-encode') {
	defined ($mpeg_encode = shift) or die $usage;
    } elsif ($arg eq '-framedir') {
	defined ($frameDir = shift) or die $usage;
    } elsif ($arg eq '-param') {
	defined ($paramFile = shift) or die $usage;
    } elsif ($arg eq '-expanded') {
	defined ($mpegFile = shift) or die $usage;
    } elsif ($arg eq '-collapsed') {
	defined ($collapsedMpegFile = shift) or die $usage;
    } elsif ($arg eq '-noexpand') {
	$buildExpanded = 0;
    } elsif ($arg eq '-nompeg') {
	$buildMpegs = 0;
    } elsif ($arg eq '-start') {
	defined ($startCol = shift) or die $usage;
    } elsif ($arg eq '-end') {
	defined ($endCol = shift) or die $usage;
    } elsif ($arg eq '-wind') {
	defined ($windPauseFlashRatio = shift) or die $usage;
    } elsif ($arg eq '-emit') {
	defined ($emitPauseFlashRatio = shift) or die $usage;
    } elsif ($arg eq '-nflash') {
	defined ($pauseTime = shift) or die $usage;
    } elsif ($arg eq '-fflash') {
	defined ($flash = shift) or die $usage;
    } elsif ($arg eq '-scale') {
	defined ($scale = shift) or die $usage;
    } elsif ($arg eq '-saveppm') {
	$savePPM = 1;
    } elsif ($arg eq '-svg') {
	$useSVG = 1;
    } elsif ($arg eq '-verbose') {
	$verbose = 1;
    } elsif ($arg =~ /^-./) {
	die $usage, "Unknown option: $arg\n";
    } else {
	push @argv, $arg;
    }
}

unless (@argv) {
    @argv = ('-');
    warn "[waiting for alignments on standard input]\n";
}

die $usage unless @argv == 1;

# do we need GD or GD::SVG?

my $imagePackage = $useSVG ? "GD::SVG" : "GD";
my $imageClass = $imagePackage . "::Image";
my $fontClass = $imagePackage . "::Font";
my $polygonClass = $imagePackage . "::Polygon";

my $imageOutputMethod = $useSVG ? "svg" : "png";
my $imageFileSuffix = $useSVG ? $svgSuffix : $pngSuffix;

eval "use $imagePackage";

# font & dimensions
my $font = $fontClass->Giant;
my $cw = $font->width;
my $ch = $font->height;

# set flash times
my ($te, $tt, $tp, $ts) = ($emitPauseFlashRatio * $pauseTime,  # emit
			   $windPauseFlashRatio * $pauseTime,  # transition (i.e. wind)
			   $pauseTime,  # pause
			   1);  # start

# read in Stockholm alignment
my $stock = Stockholm->from_file ($argv[0]);

# get tree
my $treeString = join ("", @{$stock->gf_NH});
die "Alignment must contain a phylogenetic tree\n" unless $treeString;
my $tree = Newick->parse ($treeString);

# make sure all nodes are named
for my $node (grep (!defined($tree->node_name->[$_]) || length($tree->node_name->[$_]) == 0, 0..$tree->nodes-1)) {
    my $newName = " AnonNode" . ($node + 1);  # whitespace at start of name means it won't be displayed
    while (grep (defined ($tree->node_name->[$_]) && $tree->node_name->[$_] eq $newName, 0..$tree->nodes-1)) {
	$newName .= "'";
    }
    $tree->node_name->[$node] = $newName;
    warn "Node ", $node+1, " is unnamed; calling it \"$newName\"\n";
}

# check tree topology
if ($tree->children(0) != 1) {
    # root node must be unary. Most trees won't have this, so let's help out.
    my $newRootName = " Ur-" . $tree->node_name->[0];  # whitespace at start of name means it won't be displayed
    while (grep (defined ($tree->node_name->[$_]) && $tree->node_name->[$_] eq $newRootName, 0..$tree->nodes-1)) {
	$newRootName .= "'";
    }
    $tree->insert_node_above (0, $newRootName, 0);
    warn "Added a new unary root node \"$newRootName\"\n";
}
for my $node (1..$tree->nodes - 1) {
    my $c = $tree->children($node);
    die "Internal tree node \"", $tree->node_name->[$node], "\" is not binary\n" if $c != 0 && $c != 2;
}

# get number of alignment columns
my $aligncols = $stock->columns;
die "Alignment is not flush\n" unless $stock->is_flush;

# set up startCol, endCol (Stockholm.pm uses 0-based coords; command-line args are 1-based)
$endCol = $aligncols unless defined $endCol;
--$startCol;
--$endCol;

# check all the nodes are there as sequences in the alignment
my @missing;
for (my $node = $tree->nodes - 1; $node >= 0; --$node) {
    my $node_name = $tree->node_name->[$node];
    unless (exists $stock->seqdata->{$node_name}) {
	my @child_name = map ($tree->node_name->[$_], $tree->children($node));
	die "Missing leaf sequence \"$node_name\"" if @child_name == 0;
	warn "No aligned sequence corresponding to internal tree node \"$node_name\"; estimating by principle of fewest gaps\n";
	push @{$stock->seqname}, $node_name;
	my $seqdata = "";
	for (my $col = 0; $col < $aligncols; ++$col) {
	    my $all_gaps = 1;
	    for my $child_name (@child_name) {
		if (!$stock->is_gap ($child_name, $col)) {
		    $all_gaps = 0;
		    last;
		}
	    }
	    $seqdata .= $all_gaps ? '-' : '*';
	}
	$stock->seqdata->{$node_name} = $seqdata;
    }
}

# tree topology
my @parent = @{$tree->parent};

# node depths
my @depth = (0);
for my $node (1..@parent-1) { push @depth, $depth[$parent[$node]] + 1 }
my $max_depth = max (@depth);

my @nodes_at_depth;
for my $depth (0..$max_depth) {
    push @nodes_at_depth, [grep ($depth[$_] == $depth, 0..@parent-1)];
}
my $max_nodes_at_depth = max (map (0+@$_, @nodes_at_depth));

# lay out the node positions
my @nxi = @depth;
my @nyi = map (undef, @nxi);

for (my $depth = $max_depth; $depth >= 0; --$depth) {
    my $y = 0;
    my ($last_node, $last_y);
    for my $node (@{$nodes_at_depth[$depth]}) {
	my @c = $tree->children ($node);
	if (@c == 2) {
	    $y = ($nyi[$c[0]] + $nyi[$c[1]]) / 2;  # place binary parent midway between children
	} elsif (@c == 1) {
	    $y = $nyi[$c[0]];  # place unary parent level with child
	}
	$y = max ($y, $last_y + 2) if defined($last_node) && $tree->parent->[$last_node] == $tree->parent->[$node];   # siblings need to be at least 2 apart
	$nyi[$node] = $y;
	++$y;
	$last_node = $node;
	$last_y = $y;
    }
}

# calculate total branch length of layout
sub layout_total_branch_length {
    my $len = 0;
    for my $node (1..$tree->nodes-1) {
	my $parent = $tree->parent->[$node];
	$len += sqrt (($nxi[$node]-$nxi[$parent])**2 + ($nyi[$node]-$nyi[$parent])**2);
    }
    return $len;
}

# jiggle the node positions a bit, heuristically, to try and minimize total branch length of layout
my $hard_ymax = max (@nyi);
for (my $depth = 0; $depth <= $max_depth; ++$depth) {
    my @n = @{$nodes_at_depth[$depth]};
    my @orig_nyi = @nyi;
    my @best_nyi = @nyi;
    my $best_total_len;
    my $ymin = min (map ($nyi[$_], @n));
    my $ymax = max (map ($nyi[$_], @n));
    for my $offset (-$ymin .. $hard_ymax - $ymax) {
	for my $n (@n) { $nyi[$n] = $orig_nyi[$n] + $offset }
	my $total_len = layout_total_branch_length();
	if (!defined($best_total_len) || $total_len < $best_total_len) {
	    $best_total_len = $total_len;
	    @best_nyi = @nyi;
	}
    }
    @nyi = @best_nyi;
}

# allocate width for alignment row labels
my $rowNumLen = 1 + length $tree->nodes;
my @rowName = map ($tree->node_name->[$_] . sprintf(":% ${rowNumLen}d",$_+1), 0..$tree->nodes-1);
grep (s/^\s.*:\s//, @rowName);  # remove whitespace-led names
my $labelWidth = max (map(length($_), @rowName));

# other geometric initialization
my $nVert = 12;  # node vertices

# parameters calculated by initImage()
my ($exsize, $eysize, $axsize, $aysize);
my ($xsize, $ysize, $zeroSpace);
my ($xAlignOffset, $xNodeOffset);
my (@nx, @ny);
my $nodes;
my @nodePoly;
my ($im, $white, $red, $green, $darkGreen, $blue, $cyan, $black, $purple, $midYellow);

sub initImage {
    # space for node #0
    $zeroSpace = $useZeroSpace ? $bxsize/2 + $ssize : 0;

    # EHMM dimensions
    $exsize = $bxsize * max(@nxi) + 2*$nsize + $zeroSpace + 5;
    $eysize = $bysize * max(@nyi,1) + $nsize + $ssize + $ch + 5;

    die if $exsize > 10000;

    # alignment dimensions
    $axsize = ($aligncols + $labelWidth + 1) * $cw;
    $aysize = @parent * $ch;

    # image dimensions
    $xsize = max ($exsize, $axsize);
    $ysize = $eysize + $aysize;

    # round up image to multiple of $mpegPixelMul pixels for mpeg_encode
    $xsize += $mpegPixelMul - ($xsize % $mpegPixelMul);
    $ysize += $mpegPixelMul - ($ysize % $mpegPixelMul);

    # x-offsets for alignment & tree
    $xAlignOffset = ($xsize - $axsize) / 2;
    $xNodeOffset = ($xsize - $exsize) / 2;

    # tree
    $nodes = @parent;
    @nx = map ($_*$bxsize+$nsize+$zeroSpace+$xNodeOffset, @nxi);
    @ny = map ($_*$bysize+$nsize+$ch, @nyi);

    # create a polygon for each node
    @nodePoly = ();
    foreach my $n (0..$nodes-1) {
	my $poly = $polygonClass->new;
	foreach my $s (0..$nVert-1) {
	    my $a = 2*$pi*$s/$nVert;
	    $poly->addPt ($nx[$n] + $nsize*cos($a), $ny[$n] + $nsize*sin($a));
	}
	push @nodePoly, $poly;
    }

    # allocate image and colors
    allocateImage();
}

sub allocateImage {
    $im = $imageClass->new ($xsize, $ysize);

    $white = allocateBrightRange ($im, 255, 255, 255);
    $red = allocateBrightRange ($im, 255, 0, 0);
    $green = allocateBrightRange ($im, 0, 255, 0);
    $darkGreen = allocateBrightRange ($im, 0, 127, 0);
    $blue = allocateBrightRange ($im, 0, 0, 255);
    $cyan = allocateBrightRange ($im, 0, 255, 255);
    $black = allocateBrightRange ($im, 0, 0, 0);
    $purple = allocateBrightRange ($im, 192, 0, 64);
    $midYellow = allocateBrightRange ($im, 255, 255, 128);
}

# subroutine to allocate several decreasingly bright versions of a color, and return the brightest
my $brightRange = 2;
sub allocateBrightRange {
    my ($im, $r, $g, $b) = @_;
    my @col;
    for (my $bright = $brightRange; $bright > 0; --$bright) {
	push @col, $im->colorAllocate (255 - (255-$r)*$bright/$brightRange,
				       255 - (255-$g)*$bright/$brightRange,
				       255 - (255-$b)*$bright/$brightRange);
    }
    return $col[0];
}

# flag indicating whether we're currently "recording"
my $recordFlag;

# list of PPM files generated
my @ppm;

# animation frame-drawing stuff
sub drawEHMMAlignment {
    my ($stateString, $activeNode, $alignment, $lastLen, $frame, %nodeColor) = @_;
    drawEHMM ($stateString, $activeNode, %nodeColor);
    numberNodes();
    drawAlignment ($alignment, $lastLen);
    return $recordFlag ? save ("$frameDir/$frame") : undef;
}

sub drawAlignment {
    my ($alignment, $lastLen) = @_;
    for (my $row = 0; $row < @$alignment; ++$row) {
	my $rowData = $alignment->[$row];
	my $y = $eysize + $row * $ch;
	my $rowName = $rowName[$row];
	drawString ($im, $font, $xAlignOffset + $cw * ($labelWidth - length $rowName), $y, $rowName, $black);
	if (defined $lastLen) {
	    my $lastRowLen = $lastLen->[$row];
	    my $rowLen = length($rowData);
	    if ($rowLen != $lastRowLen) {
		$im->filledRectangle ($xAlignOffset + $cw * ($labelWidth+1+$lastRowLen), $y, $xAlignOffset + $cw * ($labelWidth+1+$rowLen) - 1, $y + $ch - 1, $red)
		    if $recordFlag;
		$lastLen->[$row] = $rowLen;
	    }
	}
	drawString ($im, $font, $xAlignOffset + $cw * ($labelWidth + 1), $y, $rowData, $black);
    }
}

# subroutine to draw a string
# in PNG mode, just calls GD::Image::string once
# in SVG mode, calls it one character at a time (somehow SVG can't seem to manage fixed-width fonts...)
sub drawString {
    my ($im, $font, $x, $y, $string, $color) = @_;
    if ($recordFlag) {
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
}

# EHMM stuff
my ($estate, @align, @len, $frame, $flashCol);
my $frameListArrayRef = $buildExpanded ? [] : undef;
my $collapsedArrayRef = [];
sub initEHMM {
    initImage();
    $estate = " " . "s" x (@parent-1);
    @align = map ("", @parent);
    @len = map (0, @align);
    $frame = 0;
    mkdir $frameDir unless -e $frameDir;
    $flashCol = $black;
}

sub makeTrans {
    my ($node, $newState, $emit) = @_;
    if ($node > 1) {
	substr ($estate, $node-1, 1) = lc $newState;
    }
    if (defined $emit) {
	$align[$node-1] .= uc $emit;
    }
}

sub makeFlush {
    my $cols = max (map (length(), @align));
    foreach my $row (@align) {
	while (length ($row) < $cols) {
	    $row .= '-';
	}
    }
    @len = map ($cols, @align);
}

# drawTrans() is called by wind(), which is called by drawCol()
# calls flashFrame()
sub drawTrans {
    my ($active, $newState, $emit, $time) = @_;
    makeTrans ($active, $newState, $emit);
    flashFrame ($frameListArrayRef, $active, $time, $active-1=>$flashCol);
}

# flash on-off between a particular set of states
# flashFrame() is called by drawTrans(), drawEmit() and drawInitialFrame()
sub flashFrame {
    my ($frameListRef, $active, $time, %nodeColor) = @_;
    my $onFrame;
    if (defined $frameListRef) {
	$onFrame = drawEHMMAlignment ($estate, $active-1, \@align, \@len, $frame, %nodeColor);
	if ($recordFlag) {
	    warn "(drawing frame $frame)\n" if $verbose;
	    $frame++;
	    if (keys(%nodeColor)) {
		my $offFrame = drawEHMMAlignment ($estate, $active-1, \@align, \@len, $frame++);
		for (my $t = 0; $t < $time / 2; $t += 2) {
		    push @$frameListRef, map ($onFrame, 1..$flash), map ($offFrame, 1..$flash);
		}
	    } else {
		for (my $t = 0; $t < $time; ++$t) {
		    push @$frameListRef, $onFrame;
		}
	    }
	}
    }
    return $onFrame;
}

# frame-drawing helpers

# wind nodes to wait states
# calls drawTrans() and thus flashFrame()
sub wind {
    my @node = @_;
    foreach my $n (@node) { drawTrans ($n, "w", undef, $tt) }
}

# draw an emission, cascading down the tree
# drawEmit() is called by drawCol()
# calls flashFrame()
sub drawEmit {
    my @trans = @_;
    warn "(emitting [@trans])\n" if $verbose;
    die "Bad emit: (@trans)" if @trans % 3 != 0;
    my (@newNode, @newState, @newChar);
    for (my $i = 0; $i < @trans; $i += 3) {
	my ($node, $state, $char) = @trans[$i..$i+2];
	push @newNode, $node;
	push @newState, $state;
	push @newChar, $char;
    }
    my $emitter = $newNode[0];
    my (%nodeColor, %lastNodeColor);
    $nodeColor{$emitter-1} = $flashCol;
    makeTrans ($emitter, $newState[0], $newChar[0]);
    do {
	flashFrame ($frameListArrayRef, $emitter, $te, %nodeColor);
	my %newNodeColor;
	for (my $i = 0; $i < @newNode; ++$i) {
	    my $n = $newNode[$i] - 1;
	    my $p = $parent[$n];
	    if (exists $nodeColor{$p}) {
		$newNodeColor{$n} = $flashCol;
		makeTrans ($n+1, $newState[$i], $newChar[$i]);
	    }
	}
	%lastNodeColor = (%lastNodeColor, %nodeColor);
	%nodeColor = %newNodeColor;
    } while (keys %nodeColor);

    makeFlush();
    flashFrame ($frameListArrayRef, 0, $tp);

    flashFrame ($collapsedArrayRef, 0, $te, %lastNodeColor);
}

# draw a column, updating EHMM via wind-back and cascading emission
# calls flashFrame() via wind() and drawEmit()
sub drawCol {
    my ($column) = @_;
    my (@trans, @wind, $emitter);
    for (my $col = 0; $col < length($column); ++$col) {
	my $newChar = substr ($column, $col, 1);
	my $parentChar = $col==0 ? '-' : substr ($column, $parent[$col], 1);
	if ($newChar eq '-') {
	    if ($parentChar ne '-') {
		push @trans, $col+1, 'd', '-';
	    }
	} else {
	    $emitter = $col unless defined $emitter;
	    if ($parentChar eq '-') {
		push @trans, $col+1, 'i', $newChar;
	    } else {
		push @trans, $col+1, 'm', $newChar;
	    }
	}
    }
    for (my $n = @parent - 1; $n > $emitter; --$n) {
	my $oldState = substr ($estate, $n, 1);
	if ($oldState ne 'w') {
	    push @wind, $n+1;
	}
    }
    if (@wind) { wind (@wind) }  # calls flashFrame()
    drawEmit (@trans);  # calls flashFrame()
}

# draw initial frame
sub drawInitialFrame {
    # flashFrame() won't do anything unless given a defined array reference, so make sure it gets something even if expanded movie is disabled
    my $dummyArrayRef = $buildExpanded ? $frameListArrayRef : [];
    push @$collapsedArrayRef, flashFrame ($dummyArrayRef, 0, $ts, map (($_ => $white), 0..$nodes-1));
}

# get column text
sub getColText {
    my ($col) = @_;
    my $colText = "";
    for my $node_name (@{$tree->node_name}) {
	$colText .= substr ($stock->seqdata->{$node_name}, $col, 1);
    }
    return $colText;
}

# initialise EHMM (i.e. big transducer tree)
initEHMM();

# move EHMM to start column
$recordFlag = 0;  # don't record this
for my $col (0..$startCol-1) {
    drawCol (getColText ($col));
}

# draw initial frame
$recordFlag = 1;  # start recording
drawInitialFrame();

# loop over alignment columns, drawing frames
for my $col ($startCol..$endCol) {
    drawCol (getColText ($col));
    warn "(rendered ", $col+1-$startCol, " of ", $endCol+1-$startCol, " columns)\n" if $verbose;

    my $frameImageFile = ($frame - 1) . $imageFileSuffix;
    my $colImageFile = "$frameDir/$colImagePrefix" . ($col + 1) . $imageFileSuffix;
    symlink $frameImageFile, $colImageFile;   # make symlink to column image
}

# encode
mpegEncode ($mpegFile, $frameListArrayRef) if $buildExpanded && $buildMpegs;

# encode collapsed version
mpegEncode ($collapsedMpegFile, $collapsedArrayRef) if $buildMpegs;

# delete PPM's
warn "(deleting PPM image files)\n" if $verbose;
unless ($savePPM) {
    for my $ppm (@ppm) {
	unlink $ppm;
    }
}

# exit
exit;

# subroutine to draw states
sub drawState {
    my ($x, $y, $s, $isActive, $bgColor) = @_;  # centre of box, state type, background color
    if ($recordFlag) {  # skip all this if not recording
	my $d = $ssize*2/3;
	my $d2 = $ssize/2;
	my $dimLevel = 0;
	$im->filledRectangle ($x - $ssize, $y - $ssize, $x + $ssize, $y + $ssize, $bgColor);
	if ($s eq 'm') {  # match
	    $im->filledRectangle ($x - $ssize, $y - $d2, $x + $ssize, $y + $d2, $green + $dimLevel);
	    $im->rectangle ($x - $ssize, $y - $d2, $x + $ssize, $y + $d2, $black);
	} elsif ($s eq 'i') {  # insert
	    my $tri = $polygonClass->new;
	    $tri->addPt ($x - $ssize, $y);
	    $tri->addPt ($x + $ssize, $y + $d);
	    $tri->addPt ($x + $ssize, $y - $d);
	    $im->filledPolygon ($tri, $green + $dimLevel);
	    $im->polygon ($tri, $black);
	} elsif ($s eq 'd') {  # delete
	    my $tri = $polygonClass->new;
	    $tri->addPt ($x - $ssize, $y - $d);
	    $tri->addPt ($x - $ssize, $y + $d);
	    $tri->addPt ($x + $ssize, $y);
	    $im->filledPolygon ($tri, $darkGreen + $dimLevel);
	    $im->polygon ($tri, $black);
	} elsif ($s eq 'w') {  # wait
	    my $oct = $polygonClass->new;
	    foreach my $i (0..7) {
		my $a = ($i + .5) * $pi/4;
		$oct->addPt ($x + $d*cos($a), $y + $d*sin($a));
	    }
	    $im->filledPolygon ($oct, $red);
	    $im->polygon ($oct, $black);
	} elsif ($s eq 's') {  # start
	    my $diamond = $polygonClass->new;
	    $diamond->addPt ($x, $y + $d);
	    $diamond->addPt ($x + $d, $y);
	    $diamond->addPt ($x, $y - $d);
	    $diamond->addPt ($x - $d, $y);
	    $im->filledPolygon ($diamond, $blue + $dimLevel);
	    $im->polygon ($diamond, $black);
	} elsif ($s eq 'e') {  # end
	    $im->filledRectangle ($x - $d, $y - $d, $x + $d, $y + $d, $blue + $dimLevel);
	    $im->rectangle ($x - $d, $y - $d, $x + $d, $y + $d, $black);
	} else {
	    die "Unknown state type: $s\n";
	}
	my @b = (0..2);
	foreach my $b (@b) {
	    $im->rectangle ($x - $ssize + $b, $y - $ssize + $b, $x + $ssize - $b, $y + $ssize - $b, $black)
		if $recordFlag;
	}
	# draw a black arrow if this is the active state
	if ($isActive && $recordFlag) {
	    my ($ax, $ay) = ($x-$ssize-$anaoff, $y);
	    my $ay1;
	    if ($isActive > 0) { $ay += $ssize/2; $ay1 = $ay + $anllen }
	    else { $ay -= $ssize/2; $ay1 = $ay - $anllen }
	    arrow ($ax-$anllen, $ay1, $ax, $ay, $anlwidth, $anawidth, $analen, $black);
	}
    }
}

# filled polygon subroutine
sub autoFillPoly {
    my @point = @_;
    my $color = pop @point;
    my $poly = $polygonClass->new;
    for (my $i = 0; $i < @point; $i += 2) {
	$poly->addPt (@point[$i,$i+1]);
    }
    $im->filledPolygon ($poly, $color);
}

# subroutine to draw an arrow
sub arrow {
    my ($x1, $y1, $x2, $y2, $bwidth, $awidth, $alen, $color) = @_;
    my $dx = $x2 - $x1;
    my $dy = $y2 - $y1;
    my $norm = sqrt ($dx*$dx + $dy*$dy);
    $dx /= $norm;
    $dy /= $norm;
    my ($nx, $ny) = (-$dy, $dx);  # normal to (dx, dy)
    # back of arrow
    my $ax1 = $x2-$alen*$dx;
    my $ay1 = $y2-$alen*$dy;
    # draw line
    autoFillPoly ($x1-$bwidth*$nx, $y1-$bwidth*$ny,
		  $ax1-$bwidth*$nx, $ay1-$bwidth*$ny,
		  $ax1+$bwidth*$nx, $ay1+$bwidth*$ny,
		  $x1+$bwidth*$nx, $y1+$bwidth*$ny,
		  $color);
    # draw arrowhead
    my $ax2 = $ax1-$awidth*$nx;
    my $ay2 = $ay1-$awidth*$ny;
    my $ax3 = $ax1+$awidth*$nx;
    my $ay3 = $ay1+$awidth*$ny;
    autoFillPoly ($x2, $y2, $ax2, $ay2, $ax3, $ay3, $color);
}

# subroutine to draw branch
sub drawBranch {
    my ($x1, $y1, $x2, $y2) = @_;  # start, end of branch
    arrow ($x1, $y1, $x2, $y2, $bwidth, $awidth, $alen, $purple)
	if $recordFlag;
}

# subroutine to draw a node
sub drawNode {
    my ($n) = @_;
    $im->filledPolygon ($nodePoly[$n], $black)
	if $recordFlag;
}

# subroutine to draw a particular EHMM state
sub drawEHMM {
    my ($stateString, $activeNode, %nodeColor) = @_;
    # draw the branches and states
    foreach my $node (1..$nodes-1) {
	my $parent = $parent[$node];
	my ($nx1, $nx2) = @nx[$node,$parent];
	my ($ny1, $ny2) = @ny[$node,$parent];
	my ($sx, $sy) = (($nx1+$nx2)/2, ($ny1+$ny2)/2);
	drawBranch ($nx2, $ny2, $nx1, $ny1);
	if (length($stateString)) {
	    my $state = substr ($stateString, $node, 1);
	    drawState ($sx, $sy, $state, $node==$activeNode ? ($ny1>$ny2?1:-1) : 0, exists($nodeColor{$node}) ? $nodeColor{$node} : $white);
	}
    }
    # draw the start state
    if (length($stateString) > 0 && substr($stateString,0,1) ne " ") {
	$im->dashedLine ($ssize, $ny[0], $nx[0], $ny[0], $purple)
	    if $recordFlag;
	drawState ($ssize, $ny[0], 's', 0, $white);
    }
    # draw the nodes
    foreach my $node (0..$nodes-1) {
	drawNode ($node);
    }
}

# subroutine to clear the picture
sub clear {
    if ($useSVG) {
	allocateImage();
    } else {
	$im->filledRectangle (0, 0, $xsize, $ysize, $white);
    }
}

# subroutine to number the nodes
sub numberNodes {
    foreach my $node (0..$nodes-1) {
	drawString ($im, $font, $nx[$node] - $cw/2, $ny[$node] - $nsize - $ch, $node+1, $black);
    }
}

# subroutine to print PNG/PPM
sub save {
    my ($filename) = @_;
    my $img = "$filename$imageFileSuffix";
    my $ppm = "$filename$ppmSuffix";

    local *IMG;
    open IMG, ">$img" or die "Couldn't write to '$img': $!";
    binmode IMG;
    print IMG $im->$imageOutputMethod;
    close IMG or die "Couldn't close '$img': $!";

    system "$convert $img $ppm";
    push @ppm, $ppm;

    clear();
    return $ppm;
}

# max & min
sub max {
    my $x;
    while (@_) {
	my $y = pop @_;
	$x = $y if !defined($x) || $y > $x;
    }
    return $x;
}

sub min {
    my $x;
    while (@_) {
	my $y = pop @_;
	$x = $y if !defined($x) || $y < $x;
    }
    return $x;
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
