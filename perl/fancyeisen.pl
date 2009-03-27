#!/usr/bin/perl -w

use GD;
use Microarray::Data;

# layout options

my $block_width = 3;
my $block_height = 1;

my $name_font = gdLargeFont;
my $name_voffset = 8;

my $cluster_hspacing = 10;
my $cluster_vspacing = 10 + $name_font->height;
$cluster_vspacing = 3;

my $xcol_voffset = 2;
my $xcol_height = 4;

my $gcol_hoffset = 2;
my $gcol_width = 4;

my $bins = 128;

# parse command line

my $usage = "Usage: $0 [-bgwhite] [-columns <n>] [-blockheight <x>] [-blockwidth <x>] [-nosort] [-notext] [-allblue|-allgrey] [-rect] [-ignoresign] [-linear <min> <max>] [-rowtext] [-coltext] [-gif] <expr-file> <cluster-file> [<annot-file>]\n";

my $columns;
my $nosort = 0;
my $bgwhite = 0;
my $notext = 0;
my $allblue = 0;
my $allgrey = 0;
my $ignoresign = 0;
my $rowtext = 0;
my $coltext = 0;
my $gif = 0;
my $rectangle = 0;

my ($linmin, $linmax);
while (@ARGV && $ARGV[0] =~ /^-./) {
    my $opt = shift;
    if ($opt eq "-columns") { defined ($columns = shift) or die $usage }
    elsif ($opt eq "-blockheight") { defined ($block_height = shift) or die $usage }
    elsif ($opt eq "-blockwidth") { defined ($block_width = shift) or die $usage }
    elsif ($opt eq "-bgwhite") { $bgwhite = 1 }
    elsif ($opt eq "-nosort") { $nosort = 1 }
    elsif ($opt eq "-notext") { $notext = 1 }
    elsif ($opt eq "-allblue") { ($allblue, $allgrey) = (1, 0) }
    elsif ($opt eq "-allgrey") { ($allblue, $allgrey) = (0, 1) }
    elsif ($opt eq "-rect") { $rectangle = 1 }
    elsif ($opt eq "-ignoresign") { $ignoresign = 1 }
    elsif ($opt eq "-linear") { defined ($linmin = shift) && defined ($linmax = shift) or die $usage }
    elsif ($opt eq "-rowtext") { $rowtext = 1 }
    elsif ($opt eq "-coltext") { $coltext = 1 }
    elsif ($opt eq "-gif") { $gif = 1 }
    else { die "$usage\nUnknown option: $opt" }
}

die $usage unless @ARGV == 2 or @ARGV == 3;
my ($expr_file, $cluster_file, $annot_file) = @ARGV;

# read in files

my $dataset = Microarray::Data->new ($expr_file);

my @cluster;
my @name;
my %gene2cluster;
my %cluster_pos;
local *CLUSTER;
open CLUSTER, "<$cluster_file" or die $!;
while (<CLUSTER>) {
    my ($name, @members) = split;
    $name =~ s/_/ /g;
    push @name, $name;
    push @cluster, \@members;
    for my $i (0..@members-1) { $gene2cluster{$members[$i]} = @cluster - 1; $cluster_pos{$members[$i]} = $i }
}
close CLUSTER;

my $annot;
my %xcol;
my %gcol;
if (defined $annot_file) {
    local *ANNOT;
    open ANNOT, "<$annot_file" or die $!;
    while (<ANNOT>) {
	next unless /\S/;
	my ($directive, @arg) = split;
	if ($directive eq 'xcol') {    # directive: 'xcol <experiment-column-name> <color>'
	    my ($expt, $xcol) = @arg;
	    if (defined ($expt) && $expt ne "" && defined ($xcol) && $xcol ne "") {
		$xcol{$expt} = $xcol;
	    }
	} elsif ($directive eq 'gcol') {    # directive: 'gcol <gene-row-name> <color>'
	    my ($gene, $gcol) = @arg;
	    if (defined ($gene) && $gene ne "" && defined ($gcol) && $gcol ne "") {
		$gcol{$gene} = $gcol;
	    }
	} else {
	    warn "Unrecognised annotation directive '$directive ...'";
	}
    }
    close ANNOT;
}

# plan the layout

unless ($nosort) {
    my @i = 0..@cluster-1;
    @i = sort { @{$cluster[$a]} <=> @{$cluster[$b]} } @i;
    @cluster = @cluster[@i];
    @name = @name[@i];
}

if (!defined $columns) { $columns = int sqrt (@cluster + 0) }

my $max_rowtext_len = 0;
my $max_coltext_len = 0;
foreach my $row_name (@{$dataset->row_name}) {
    if (length($row_name) > $max_rowtext_len) { $max_rowtext_len = length($row_name) }
}
foreach my $col_name (@{$dataset->column_name}) {
    if (length($col_name) > $max_coltext_len) { $max_coltext_len = length($col_name) }
}
my $rowtext_offset = $rowtext ? $name_font->width * ($max_rowtext_len + .5) : 0;
my $coltext_offset = $coltext ? $name_font->height * 1.5 : 0;

my @row_y_orig;
my @row_vblocks;
my $img_height = 0;
for (my $c = 0; $c < @cluster; $c += $columns) {
    my $max_blocks;
    for (my $i = 0; $i < $columns && $i + $c < @cluster; ++$i) {
	my $blocks = @{$cluster[$i + $c]};
	if (!defined ($max_blocks) || $blocks > $max_blocks) { $max_blocks = $blocks }
    }
    $img_height += $cluster_vspacing + $coltext_offset;
    push @row_y_orig, $img_height;
    push @row_vblocks, $max_blocks;
    $img_height += $block_height * $max_blocks;
}
my $rows = @row_y_orig + 0;

my $img_width = $columns * $block_width * $dataset->columns + $columns * ($cluster_hspacing + $rowtext_offset);

# create the image

my $img = GD::Image->new ($img_width, $img_height);

# generate red <--> green color scale

my ($black, $white);
sub make_black { shift->colorAllocate (0, 0, 0) }
sub make_white { shift->colorAllocate (255, 255, 255) }
if ($bgwhite) { $white = make_white($img); $black = make_black($img) }
else { $black = make_black($img); $white = make_white($img) }

my $font_col = $bgwhite ? $black : $white;

my $dim = $bgwhite ? 255-16 : 16;
my $grey = $img->colorAllocate ($dim, $dim, $dim);

my ($yellow, $green, $blue, $red, $orange);
my %annot_color = ( 'y' => $yellow = $img->colorAllocate (255, 255, 0),
		    'g' => $green = $img->colorAllocate (0, 255, 0),
		    'b' => $blue = $img->colorAllocate (0, 0, 255),
		    'r' => $red = $img->colorAllocate (255, 0, 0),
		    'o' => $orange = $img->colorAllocate (255, 127, 0) );

sub make_red {
    my ($img, $x) = @_;
    $bgwhite ? $img->colorAllocate (255, 255-$x, 255-$x) : $img->colorAllocate ($x, 0, 0);
}

sub make_green {
    my ($img, $x) = @_;
    $bgwhite ? $img->colorAllocate (255-$x, 255, 255-$x) : $img->colorAllocate (0, $x, 0);
}

sub make_blue {
    my ($img, $x) = @_;
    $bgwhite ? $img->colorAllocate (255-$x, 255-$x, 255) : $img->colorAllocate (0, 0, $x);
}

sub make_grey {
    my ($img, $x) = @_;
    $img->colorAllocate (255*(1-$x), 255*(1-$x), 255*(1-$x));
}

my @color;
if ($allgrey) {
    for (my $i = 0; $i < $bins; ++$i) { push @color, make_grey ($img, $i / $bins) }
} elsif ($allblue) {
    for (my $i = 0; $i < $bins; ++$i) { push @color, make_blue ($img, $i) }
} else {
    for (my $i = int($bins/2); $i > 0; --$i) { push @color, make_red ($img, 2*$i) }
    for (my $i = 0; $i < int($bins/2); ++$i) { push @color, make_green ($img, 2*$i) }
}

if ($dataset->min >= 0) {
    $dataset->offset (($dataset->min-$dataset->max)/2);
}
my @bin;
if (defined $linmin) { @bin = map ($linmin + $_*($linmax-$linmin)/($bins-1), 0..$bins-1) }
elsif ($ignoresign) { @bin = $dataset->get_bins_ignore_sign ($bins) }
else { @bin = $dataset->get_bins ($bins) }

# subroutine to find the right color bin by binary chop

sub get_bin {
    my ($bin_ref, $val) = @_;
    my $min_b = 0;
    my $step_b = @$bin_ref - 1;
    while ($step_b > 0) {
	my $test_b = $min_b + $step_b;
	if ($test_b < @$bin_ref && $$bin_ref[$test_b] < $val) {
	    $min_b = $test_b;
	} else {
	    $step_b = int ($step_b / 2);
	}
    }
    return $min_b;
}

# subroutine to draw a row

sub draw_row {
    my ($x_orig, $y_orig, $row_data, $row_name, $dataset) = @_;
    my ($min, $max) = ($dataset->min, $dataset->max);
    for (my $column = 0; $column < @$row_data; ++$column) {
	my $x = $x_orig + $block_width * $column;
	my $data = $$row_data[$column];
	my $color = $grey;
	if (defined($data) && $data =~ /\d/) {
	    my $color_index = get_bin (\@bin, $data);
	    $color = $color [$color_index];
	    
	    if ($rectangle) {
		$img->filledRectangle ($x,
				       $y_orig,
				       $x + $block_width - 1,
				       $y_orig + $block_height - 1,
				       $color);
	    } else {
		my $cx = $x + $block_width / 2;
		my $cy = $y_orig + $block_height / 2;

		my $scaled = ($data - $min) / ($max - $min);
		my $radius = sqrt ($scaled);

		my $rw = $radius * $block_width;
		my $rh = $radius * $block_height;
		$rw = 1 if $rw < 1;
		$rh = 1 if $rh < 1;

		$color = $color[@color-1];    # hack

		$img->arc ($cx, $cy, $rw, $rh, 0, 360, $black);
		$img->fill ($cx, $cy, $color);
	    }
	} else {
#	    $img->line ($x+1,
#			$y_orig+1,
#			$x + $block_width - 2,
#			$y_orig + $block_height - 2,
#			$black);
#	    $img->line ($x + $block_width - 2,
#			$y_orig+1,
#			$x+1,
#			$y_orig + $block_height - 2,
#			$black);
	    $img->line ($x + 1,
			$y_orig + $block_height / 2,
			$x + $block_width - 2,
			$y_orig + $block_height / 2,
			$black);
	}
    }
    # draw gene label
    if ($rowtext) {
	$img->string ($name_font, $x_orig - $name_font->width * (length($row_name) + .5), $y_orig, $row_name, $font_col);
    }
}

# draw the clusters

for (my $row = 0; $row < $rows; ++$row) {
    for (my $column = 0; $column < $columns; ++$column) {
	my $x_orig = $column * $block_width * $dataset->columns + ($column + 1) * $cluster_hspacing + $rowtext_offset;
	my $c = $row * $columns + $column;
	last if $c >= @cluster;
	my $cref = $cluster[$c];
	my $y_top = $row_y_orig[$row];
	for (my $i = 0; $i < @$cref; ++$i) {
	    my $y_orig = $y_top + $i * $block_height;
	    draw_row ($x_orig, $y_orig, $dataset->get_row_by_name ($$cref[$i]), $$cref[$i], $dataset);
	}
	# draw colored annotation bars for each experiment
	while (my ($expt, $xcol) = each %xcol) {
	    my $x = $x_orig + $block_width * $dataset->column_index->{$expt};
	    my $color = $annot_color{$xcol};
	    $img->filledRectangle ($x,
				   $y_top - $xcol_voffset - $xcol_height + 1,
				   $x + $block_width - 1,
				   $y_top - $xcol_voffset,
				   $color);
	}
	# draw experiment label
	for (my $e = 0; $e < @{$dataset->column_name}; ++$e) {
	    my $expt = $dataset->column_name->[$e];
	    my $x = $x_orig + $block_width * ($e + .5) - $name_font->width * length($expt) / 2;
	    my $y = $y_top - $xcol_voffset - $xcol_height - $name_font->height;
	    $img->string ($name_font, $x, $y, $expt, $font_col);
	}
	# draw colored annotation bars for each gene
	while (my ($gene, $gcol) = each %gcol) {
	    next unless $gene2cluster{$gene} == $c;
	    my $y = $y_top + $cluster_pos{$gene} * $block_height;
	    my $color = $annot_color{$gcol};
	    $img->filledRectangle ($x_orig - $gcol_hoffset - $gcol_width + 1,
				   $y,
				   $x_orig - $gcol_hoffset,
				   $y + $block_height - 1,
				   $color);
	}
	# draw cluster name
	unless ($notext) {
	    my $name = $name[$c];
	    my $x = $x_orig + $block_width * $dataset->columns / 2;
	    $x -= $name_font->width * length($name) / 2;
	    my $y = $y_top - $name_voffset - $name_font->height * 2.5;
	    $img->string ($name_font, $x, $y, $name, $font_col);
	}
    }
}

# output

binmode STDOUT;
if ($gif) { print $img->gif }
else { print $img->png }
