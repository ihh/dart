#!/usr/bin/env perl -w

my $cellsize = 7;
my $border = 2;
my $font = "Courier";

local $_;

$_ = <>;
my ($c, $a) = split;

$_ = <>;
my @alph = split;

$_ = <>;
my @eqm = split;

my @rate;
my $max = 0;
for (my $i = 0; $i < @alph; ++$i) {
    $_ = <>;
    my @row = split;
    die "Row $i has ", @row+0, " elements, not ", @alph+0 unless @row == @alph;
    push @rate, \@row;
    for (my $j = 0; $j < @alph; ++$j) {
	$max = $row[$j] if $i != $j && $row[$j] > $max;
    }
}

my $a1 = $a - 1;
my $maxsymlen = 0;
foreach my $sym (@alph) {
    $maxsymlen = length ($sym) if length ($sym) > $maxsymlen;
}

my $labsize = $cellsize * $maxsymlen;
my $xoffset = $labsize + $border + 1;
my $yoffset = $border + 1;

my $xsize = $xoffset + $cellsize * $a + $border;
my $ysize = $yoffset + $cellsize * $a + $border + $labsize;

$xsize += 20; $ysize += 20;   # hack!

print map ("$_\n",
	   "\%!PS-Adobe-2.0",
	   "\%\%BoundingBox: 0 0 $xsize $ysize",
	   "/alphsize $a def",
	   "/border $border def",
	   "/cellsize $cellsize def",
	   "/font ($font) def",
	   "/dot {",
	   "/area exch def",
	   "/y exch alphsize 1 sub exch sub def",
	   "/x exch def",
	   "newpath",
	   "x 0.5 add cellsize 1 add mul y 0.5 add cellsize 1 add mul",
	   "area sqrt cellsize mul 2 div 0 360 arc fill",
	   "} def",
	   "/xlabel {",
	   "/text exch def",
	   "/y exch alphsize 1 sub exch sub def",
	   "/x 0 def",
	   "gsave",
	   "font findfont cellsize scalefont setfont",
	   "x cellsize 1 add mul y cellsize 1 add mul moveto",
	   "text stringwidth pop border add 2 add -1 mul 0 rmoveto",
	   "text show",
	   "grestore",
	   "} def",
	   "/ylabel {",
	   "/text exch def",
	   "/y $a def",
	   "/x exch def",
	   "gsave",
	   "font findfont cellsize scalefont setfont",
	   "x cellsize 1 add mul y cellsize 1 add mul moveto",
	   "90 rotate",
	   "border 2 add cellsize -1 mul rmoveto",
	   "text show",
	   "grestore",
	   "} def",
	   "$xoffset $yoffset translate",
	   "gsave",
	   "$border setlinewidth 0.5 setgray",
	   "newpath",
	   "border -2 div -1 add dup moveto",
	   "$a cellsize 1 add mul border 2 add add dup dup dup 0 rlineto",
	   "0 exch rlineto",
	   "-1 mul 0 rlineto",
	   "-1 mul 0 exch rlineto",
	   "closepath stroke",
	   "grestore");

for (my $i = 0; $i < @alph; ++$i) {
    my $alph = uc $alph[$i];
    print "$i ($alph) xlabel $i ($alph) ylabel\n";
    for (my $j = 0; $j < @alph; ++$j) {
	if ($i != $j) {
	    printf "%d %d %.3g dot\n", $j, $i, $rate[$i]->[$j] / $max;
	}
    }
}

print "showpage\n", "\%\%EOF\n";
