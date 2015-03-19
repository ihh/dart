#!/usr/bin/env perl -w

# set up pagers
my $less = "less -r";
my $cat = "cat";
my $pipe_to_cat = 1;

# Stockholm tags & special characters
my $PS = "PS";  # the #=GC tag for primary sequence
my $SS = "SS";  # the #=GC tag for secondary structure
my $SC = "SC";  # the #=GC tag for score
my $PS_CONS = "PS_cons";  # the #=GC tag for primary sequence consensus
my $SS_CONS = "SS_cons";  # the #=GC tag for secondary structure consensus
my $SN = "SN";  # the #=GC tag for stem number
my $BP = "BP";  # the #=GC tag for number of basepaired rows
my $CS = "CS";  # the #=GC and #=GF tag for number of covariant substitutions

my %rchar = ('('=>')','<'=>'>','['=>']');
my %lchar = map (($rchar{$_}=>$_), keys(%rchar));
my %gapchar = map (($_=>1), '-', '.', '_', ',');

my %strongly_bonded = map (($_ => 1), qw(AU GC CG UA GU UG));

# display parameters
my $screenColumns = (`tput cols` + 0) || 80;
my $covOnly = 0;
my $canOnly = 0;
my $raw = 0;
my $diffName;
my $refName;

# stem characters
my @stemChar = (0..9, map (chr($_), ord('A')..ord('Z'), ord('a')..ord('Z')));
my $stemChars = @stemChar;

# command-line args
my $summOnly = 0;
my $html = 0;
my %invert;

# help message
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname <Stockholm file>\n";
$usage .=   "          [-h] print this message\n";
$usage .=   "       [-less] pipe file through '$less' instead of '$cat'\n";
$usage .=   "   [-cols <N>] number of columns per row (default is $screenColumns)\n";
$usage .=   "   [-gc <TAG>] '#=GC' tag for secondary structure; default is '$SS_CONS'\n";
$usage .=   "        [-cov] only color columns with compensatory mutations\n";
$usage .=   "        [-can] only color canonical (and wobble) base-pairs\n";
$usage .=   "  [-inv <SEQ>] invert background/foreground for SEQ\n";
$usage .=   " [-diff <SEQ>] invert background/foreground for mutations relative to SEQ\n";
$usage .=   "  [-ref <SEQ>] use alternate color scheme, with SEQ as reference (see perldoc)\n";
$usage .=   "       [-summ] don't print alignment, just print covariance summary\n";
$usage .=   "                ('covariance' here means 'having compensatory mutations')\n";
$usage .=   "        [-raw] don't print $CS or $BP annotation lines\n";
$usage .=   "       [-html] print results in HTML\n";
$usage .=   "\n";
$usage .=   "For more help, type the following:\n";
$usage .=   "perldoc $0\n";
$usage .=   "\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-less") {
	$pipe_to_cat = 0;
    } elsif ($arg eq "-cols") {
	defined ($screenColumns = shift) or die $usage;
    } elsif ($arg eq "-gc") {
	defined ($SS_CONS = shift) or die $usage;
    } elsif ($arg eq "-cov") {
	$covOnly = 1;
    } elsif ($arg eq "-can") {
	$canOnly = 1;
    } elsif ($arg eq "-diff") {
	defined ($diffName = shift) or die $usage;
    } elsif ($arg eq "-inv") {
	my $seqname;
	defined ($seqname = shift) or die $usage;
	$invert{$seqname} = 1;
    } elsif ($arg eq "-ref") {
	defined ($refName = shift) or die $usage;
    } elsif ($arg eq "-summ") {
	$summOnly = 1;
    } elsif ($arg eq "-raw") {
	$raw = 1;
    } elsif ($arg eq "-html") {
	$html = 1;
    } else {
	push @argv, $arg;
    }
}
unless (@argv) {
    @argv = ('-');
    warn "[waiting for alignments on standard input]\n";
}
die $usage unless @argv == 1;

# check for conflicting options
if (defined($refName) && ($covOnly || $canOnly || defined($diffName) || keys(%invert) > 0)) {
    die "\nERROR: -ref option cannot be used with -cov, -can, -diff or -inv\n\n";
}

# get input filename
my ($stockfile) = @argv;

# set up the color codes
my ($header, $footer, $white, $invWhite, @ansi, @invAnsi);
if ($html) {
# Define the document style sheet
    $header = make_header();
    $footer = "</span></pre></body></html>\n";

    $white = '</span><span class="W">';
    $invWhite = '</span><span class="IW">';
    @ansi = map ("</span><span class=\"$_\">", qw(R G Y B M C IR IG IY IB IM IC IW));
    @invAnsi = map ("</span><span class=\"$_\">", qw(IR IG IY IB IM IC R G Y B M C W));
} else {
    $header = $footer = "";
    $white = e(27).e(37);  # white text
    $invWhite = e(7).e(37);  # white text
    @ansi = (map (e(27).e(30+$_), 1..6),  # non-inverse, red..cyan
	     map (e(7).e(30+$_), 1..7));  # inverse, red..white
    @invAnsi = (map (e(7).e(30+$_), 1..6),  # inverse, red..cyan
		map (e(27).e(30+$_), 1..7));  # non-inverse, red..white
}

# alternate, sscolor-like color scheme
my %altColor = ('noPair' => -1,  # white
		'noPairSub' => 4,  # magenta
		'compatSgl' => 3,  # blue
		'compatDbl' => 1,  # green
		'nonCompat' => 0,  # red
		'noSub' => -1);  # white

# HTML



# other color-related stuff
my $maxColor = @ansi;


# more Stockholm stuff
my $ps_cons_hdr = "#=GC $PS_CONS";
my $ss_cons_hdr = "#=GC $SS_CONS";
my $sn_hdr = "#=GC $SN";
my $sc_hdr = "#=GF $SC";
my $cs_hdr = "#=GF $CS";
my $gc_bp_hdr = "#=GC $BP";
my $gc_cs_hdr = "#=GC $CS";

# open less pager
my $pager = $pipe_to_cat ? $cat : $less;
local *LESS;
open (LESS, "|$pager") or open (LESS, ">-");

# print header
print LESS $header;

# OK, let's go
local *STOCK;
open STOCK, "<$stockfile" or die "Couldn't open '$stockfile': $!";

while (1) {
    my (@seqname, %seq, %ss, %gc, %gr, $score);
    my ($ps_cons, $ss_cons, $preamble);
    my ($basepairs, $stems, $covStems) = (0, 0, 0);
    my $seqlen = 0;
    my $found_separator = 0;
    while (<STOCK>) {
	if (/^\s*\#=GR\s*(\S+)\s*$SS\s*(\S+)\s*$/) { $ss{$1} = "" unless exists $ss{$1}; $ss{$1} .= $2 }  # by-seq secondary structure
	elsif (/^\s*\#=GR\s*(\S+)\s*(\S+)\s*(\S+)\s*$/) { $gr{$2} = {} unless exists $gr{$1}; $gr{$1}->{$2} = "" unless exists $gr{$1}->{$2}; $gr{$1}->{$2} .= $3 }  # miscellaneous #=GR lines
	elsif (/^\s*\#=GC\s*$PS_CONS\s*(\S+)\s*$/) { $ps_cons = "" unless defined $ps_cons; $ps_cons .= $1 }  # consensus primary sequence
	elsif (/^\s*\#=GC\s*$SS_CONS\s*(\S+)\s*$/) { $ss_cons = "" unless defined $ss_cons; $ss_cons .= $1 }  # consensus secondary structure
	elsif (/^\s*\#=GC\s*(\S+)\s*(\S+)\s*$/) { $gc{$1} = "" unless exists $gc{$1}; $gc{$1} .= $2 }  # miscellaneous #=GC lines
	elsif (/^\s*\#=GF\s*$SC\s*(\S+)\s*$/) { $score = $1 }  # score
	elsif (/^\s*\#/) { $preamble = "" unless defined $preamble; $preamble .= $_ }  # unrecognised line starting with '#'; append to preamble
	elsif (/^\s*\/\//) { $found_separator = 1; last }  # alignment separator
	elsif (/^\s*(\S+)\s*(\S+)\s*$/) { unless (exists $seq{$1}) { $seq{$1} = ""; push @seqname, $1 } $seq{$1} .= $2; $seqlen = length($seq{$1}) }
	elsif (/\S/) { warn "Ignoring line: $_" }
    }
    if (@seqname) {
	my @stem = map (-1, 1..$seqlen);  # stem index for each column
	my @cov = map (0, 1..$seqlen);  # flag indicating whether column is in a covariant basepair
	my %canByRow = map (($_ => [@cov]), @seqname);  # hash indicating whether each row-column position is in a canonical basepair
	my %subByRow;  # hash recording number of substitutions relative to reference sequence
	if (defined $refName) {
	    if (!exists $seq{$refName}) {
		warn "WARNING: reference sequence '$refName' not found in alignment\n";
		$seq{$refName} = "N" x $seqlen;
	    }
	    my $refSeq = $seq{$refName};
	    while (my ($seqname, $seq) = each %seq) { $subByRow{$seqname} = [map (uc(substr($refSeq,$_,1)) eq uc(substr($seq,$_,1)) ? 0 : 1, 0..length($seq)-1)] }
	}
	if (defined $diffName) {
	    if (!exists $seq{$refName}) {
		warn "WARNING: reference sequence '$diffName' not found in alignment\n";
		$seq{$diffName} = "N" x $seqlen;
	    }
	}
	my @bp = map (0, 1..$seqlen);  # number of actual basepaired rows in each basepaired column-pair
	my @nbr;  # neighbours of each stem
	if (defined($ss_cons) && length($ss_cons)) {
	    # figure out colors from SS_CONS row
	    my $currentStem = -1;
	    my ($lastl, $lastr);
	    my $stemcol = 0;
	    my @lpos;
	    for (my $i = 0; $i < length($ss_cons); ++$i) {
		# denote this stem by '<>', neighbours by '()' and '[]'
		my $c = substr ($ss_cons, $i, 1);
		if (exists $rchar{$c}) {  # <
		    push @lpos, $i;
		} elsif (exists $lchar{$c}) {  # >
		    unless (@lpos) { die "Too many >'s in line: $ss_cons_hdr $ss_cons\n" }
		    my $l = pop @lpos;
		    my $r = $i;
		    ++$basepairs;
		    # is column covariant?
		    my %pairs;
		    my $bpTot = 0;
		    while (my ($seqname, $seq) = each %seq) {
			my $lc = uc substr ($seq, $l, 1);
			my $rc = uc substr ($seq, $r, 1);
			my $pair = $lc.$rc;
			$pair =~ s/T/U/g;
			if (exists $strongly_bonded{$pair}) {
			    ++$pairs{$pair};
			    ++$bpTot;
			    $canByRow{$seqname}->[$l] = 1;
			    $canByRow{$seqname}->[$r] = 1;
			    if (defined $refName) {
				my $lcRef = uc substr ($seq{$refName}, $l, 1);
				my $rcRef = uc substr ($seq{$refName}, $r, 1);
				# give the RIGHT subByRow count an extra nudge if the LEFT char doesn't match, and vice versa
				++$subByRow{$seqname}->[$l] if $rc ne $rcRef;
				++$subByRow{$seqname}->[$r] if $lc ne $lcRef;
			    }
			}
		    }
		    my $nCov = keys(%pairs);
		    if ($nCov > 1) {
			$cov[$l] = $cov[$r] = $nCov;
			++$covStems;
		    }
		    $bp[$l] = $bp[$r] = $bpTot;
		    # figure out stem numbering
		    my $first = !defined($lastl) || !defined($lastr);
		    if ($first || $lastl != $l+1 || $lastr != $r-1) {  # excluding  <()>   possibilities are   <().>   <.()>   <.().>   ()<>
			++$currentStem;
			my @currentNbr;
			if (!$first) {
			    my $prevStem = $stem[$lastl];
			    push @currentNbr, $prevStem;
			    push @{$nbr[$prevStem]}, $currentStem;
			    # catch the possibility of  <[]()>
			    for (my $j = $l+1; $j < $lastl; ++$j) {
				my $jStem = $stem[$j];
				if ($jStem >= 0 && $jStem != $prevStem) {
				    push @{$nbr[$jStem]}, $currentStem;
				    push @currentNbr, $jStem;
				    last;
				}
			    }
			}
			push @nbr, \@currentNbr;
			++$stems;
		    }
		    $stem[$l] = $stem[$r] = $currentStem;
		    ($lastl, $lastr) = ($l, $r);
		}
	    }
	    if (@lpos) { die "Too many <'s in $SS_CONS string: $ss_cons\n" }
	}

	# make summary line
	my $summary = "$seqlen columns, $basepairs paired, $covStems covariant; $stems stems";

	# either print summary line, or proceed to coloring & display
	if ($summOnly) {
	    print LESS $summary, "\n";
	} else {

	    # figure out stem colors so neighbours always different
	    my @stemColor = map (undef, @nbr);  # color of each stem
	    my @stemNum = map (undef, @nbr);  # recode stem numbers in sensible order
	    my $currentColor = -1;
	    my $currentStemNum = -1;
	    my $lasts = -1;
	    for (my $i = 0; $i < $seqlen; ++$i) {
		my $s = $stem[$i];
		if ($s >= 0 && $s != $lasts) {
		    if (!defined $stemColor[$s]) {
			while (1) {
			    $currentColor = ($currentColor + 1) % $maxColor;
			    ++$currentStemNum;
			    my $unique = 1;
			    foreach my $nbr (@{$nbr[$s]}) {
				if (defined($stemColor[$nbr]) && $stemColor[$nbr] == $currentColor) {
				    $unique = 0;
				    last;
				}
			    }
			    last if $unique;
			}
			$stemColor[$s] = $currentColor;
			$stemNum[$s] = $currentStemNum;
		    }
		}
		$lasts = $s;
	    }

	    # figure out name width
	    my (@gr_hdr, @ss_hdr);
	    while (my ($name, $hash) = each %gr) { push @gr_hdr, map ("#=GR $name $_", keys %$hash) }
	    for my $name (@seqname) {
		push @ss_hdr, exists($ss{$name}) ? "#=GR $name $SS" : $name;
	    }
	    my $w = max (map (length, $ss_cons_hdr, @ss_hdr, @gr_hdr));
	    while (my ($name, $gr) = each %gr) {
		$w = max ($w, map (length, map ("#=GR $name $_", keys %$gr)));
	    }

	    # figure out display columns
	    my $seqColumns = $screenColumns - $w - 1;

	    # number the stems
	    my $stemString = join ("", map ($_<0?'.':$stemChar[$stemNum[$_] % $stemChars], @stem));

	    # flag covariant columns
	    my $covChars = join ("", map ($_ ? sprintf("%x",$_) : ".", @cov));

	    # count basepairs in each column
	    my @bpCharMap = ('.', 1..9, map(chr($_), ord('a')..ord('z')), '+');
	    my $bpChars = join ("", map ($_ >= @bpCharMap ? $bpCharMap[@bpCharMap-1] : $bpCharMap[$_], @bp));

	    # hack the @stemColor array if we're using the alternate color scheme
	    my $stemColRef = defined($refName) ? [0..7] : \@stemColor;

	    # print the alignment
	    print LESS $preamble if defined $preamble;
	    for (my $col = 0; $col < $seqlen; $col += $seqColumns) {

		my $width = min ($seqColumns, $seqlen - $col);
		my @printRowArgs = ($col, $width, $w);

		my @stemSlice = @stem[$col .. $col + $width - 1];

		my (@colAnsi, @covColAnsi);
		if (defined $refName) {

		    @covColAnsi = @colAnsi = map ($white, @stemSlice);

		} else {

		    @colAnsi = map (getStemColor($_,$stemColRef), @stemSlice);

		    @covColAnsi =
			$covOnly
			? map ($cov[$_+$col] ? $colAnsi[$_] : $white, 0..@colAnsi-1)
			: @colAnsi;

		}

		printRow ($white, $gc_cs_hdr, $covChars, @printRowArgs, \@covColAnsi) unless $raw;
		printRow ($white, $gc_bp_hdr, $bpChars, @printRowArgs, \@covColAnsi) unless $raw;
		for (my $i = 0; $i < @seqname; ++$i) {

		    my $name = $seqname[$i];
		    my $invertFlag = exists($invert{$name}) ? 1 : 0;

		    my $nameColor = $invertFlag ? $invWhite : $white;

		    my @rowStemSlice;
		    if (defined $refName) {

			@rowStemSlice = map ($stemSlice[$_] < 0  # single-stranded?
					     ? ($subByRow{$name}->[$col+$_]  # any substitutions?
						? $altColor{'noPairSub'}  # single-stranded, one substitution
						: $altColor{'noPair'})  # single-stranded, no substitution
					     : ($canByRow{$name}->[$col+$_]  # canonical?
						? ($subByRow{$name}->[$col+$_] == 2 ? $altColor{'compatDbl'}  # base-paired, canonical, two substitutions
						   : ($subByRow{$name}->[$col+$_] == 1 ? $altColor{'compatSgl'}  # base-paired, canonical, one substitution
						      : $altColor{'noSub'}))  # base-paired, canonical, no substitutions
						: ($subByRow{$name}->[$col+$_]  # any substitutions?
						   ? $altColor{'nonCompat'}  # base-paired, non-canonical, one or more substitutions
						   : $altColor{'noSub'})),  # base-paired, non-canonical, no substitutions
					     0..@stemSlice-1);

		    } else {
			@rowStemSlice = @stemSlice;
			@rowStemSlice = map ($cov[$col+$_] ? $rowStemSlice[$_] : -1, 0..@rowStemSlice-1) if $covOnly;
			@rowStemSlice = map ($canByRow{$name}->[$col+$_] ? $rowStemSlice[$_] : -1, 0..@rowStemSlice-1) if $canOnly;
		    }

		    my @invertFlag = map ($invertFlag, @rowStemSlice);
		    @invertFlag = map (uc(substr($seq{$name},$col+$_,1)) eq uc(substr($seq{$diffName},$col+$_,1)) ? $invertFlag[$_] : !$invertFlag[$_], 0..@invertFlag-1) if defined($diffName) && $diffName ne $name;

		    my @rowColAnsi = map (getStemColor($_,$stemColRef), @rowStemSlice);

		    my @invColAnsi = map ($invertFlag[$_] ? getInverseStemColor($rowStemSlice[$_],$stemColRef) : getStemColor($rowStemSlice[$_],$stemColRef), 0..@rowStemSlice-1);

		    printRow ($nameColor, $name, $seq{$name}, @printRowArgs, \@invColAnsi);
		    printRow ($white, $ss_hdr[$i], $ss{$name}, @printRowArgs, \@covColAnsi) if exists $ss{$name};
		    if (exists $gr{$name}) {
			foreach my $tag (sort keys %{$gr{$name}}) {
			    printRow ($white, "#=GR $name $tag", $gr{$name}->{$tag}, @printRowArgs, \@covColAnsi);
			}
		    }
		}
		printRow ($white, $ps_cons_hdr, $ps_cons, @printRowArgs, \@colAnsi) if defined $ps_cons;
		printRow ($white, $ss_cons_hdr, $ss_cons, @printRowArgs, \@colAnsi) if defined $ss_cons;
		foreach my $tag (sort keys %gc) {
		    printRow ($white, "#=GC $tag", $gc{$tag}, @printRowArgs, \@colAnsi);
		}
		printRow ($white, $sn_hdr, $stemString, @printRowArgs, \@colAnsi) if defined($ss_cons) && length($ss_cons);
		print LESS "\n";
	    }
	    printRowBase ($white, $sc_hdr, $score, $w) if defined $score;
	    printRowBase ($white, $cs_hdr, $summary, $w) unless $raw;
	    print LESS "//\n";
	}
    }
    # if no alignment separator, then quit the loop
    last unless $found_separator;
}
close STOCK;

# print footer
print LESS $footer;
close LESS;

# max subroutine
sub max {
    my ($x, @y) = @_;
    foreach my $y (@y) { $x = $y if $y > $x }
    return $x;
}

# min subroutine
sub min {
    my ($x, @y) = @_;
    foreach my $y (@y) { $x = $y if $y < $x }
    return $x;
}

# escape code subroutine
sub e {
    my $code = shift;
    return chr(27)."[$code"."m";
}

# column-coloring subroutine.
# @$colAnsi is an array with length($string) elements,
# each element corresponding to a string of ANSI terminal characters to go before the corresponding subchar of $string
sub colorString {
    my ($string, $colAnsi) = @_;
    my $len = length $string;
    return $string if $len != @$colAnsi;
    my $out = "";
    my $lastCol = "";
    for (my $i = 0; $i < $len; ++$i) {
	my $nextCol = $$colAnsi[$i];
	if ($nextCol ne $lastCol) { $out .= $lastCol . $nextCol; $lastCol = $nextCol }
	$out .= substr ($string, $i, 1);
    }
    $out .= $lastCol;
    return $out;
}

# subroutines to get color for a particular stem
sub getStemColor {
    my ($stemNum, $stemColorRef) = @_;
    return $stemNum < 0 ? $white : $ansi[$stemColorRef->[$stemNum]];
}

sub getInverseStemColor {
    my ($stemNum, $stemColorRef) = @_;
    return $stemNum < 0 ? $invWhite : $invAnsi[$stemColorRef->[$stemNum]];
}

# base row printing subroutine
sub printRowBase {
    my ($hdrColor, $hdr, $row, $nameWidth) = @_;
    print LESS "$hdrColor$hdr$white", " " x ($nameWidth + 1 - length($hdr)), $row, $white, "\n";
}

# colored row printing subroutine
# @$colAnsi is an array with length($row) elements,
# each element corresponding to a string of ANSI terminal characters to go before the corresponding subchar of $row
sub printRow {
    my ($hdrColor, $hdr, $row, $startCol, $nCols, $nameWidth, $colAnsi) = @_;
    my $rowlen = length $row;
    if ($startCol >= $rowlen) {
	$startCol = 0;
	$nCols = 0;
    } else {
	$nCols = min ($nCols, $rowlen - $startCol);
    }
    my $rowSubstr = substr ($row, $startCol, $nCols);
    printRowBase ($hdrColor, $hdr, colorString ($rowSubstr, $colAnsi), $nameWidth);
}

=head1 NAME

colorstock.pl - colorize RNA multiple alignments in ANSI terminals

=head1 SYNOPSIS

 colorstock.pl -h

 colorstock.pl FILENAME
               [-less] [-cols N] [-gc TAG] [-cov] [-can]
               [-inv SEQ] [-ref SEQ] [-diff SEQ]
               [-summ] [-html] 

 cat FILENAME1 [FILENAME2 ...] | colorstock.pl [options]

=head1 DESCRIPTION

The program uses ANSI terminal color (or, optionally, HTML) to identify unbroken stems in a multiple alignment of noncoding RNA (or in a series of such alignments)
and to highlight compensatory mutations or individual rows in the alignment.

The input alignments must be in Stockholm format.
The output contains nonprintable control characters that will be interpreted as colors by an ANSI-compliant terminal (including e.g. VT100 or X11 terminals).
If the HTML output is selected, then colors can be viewed using a standards-compliant web browser.

For the purposes of this program, an "unbroken stem" is defined to be a stacked sequence of base-pairs uninterrupted by interior loops or bulges.
Any such interruptions will result in a new stem being counted and a new color being used for the resumed stem.

A column is defined to display a compensatory mutation if (i) it is base-paired with another column,
(ii) when considering all pairs of nucleotides in all rows at those two column positions, there are at least two different Watson-Crick or wobble base-pairs.

=head1 INPUT FORMAT

The input file must be in Stockholm format, defined here:

L<http://biowiki.org/StockholmFormat>

A list of utilities for format conversions and other operations on Stockholm files can be found here:

L<http://biowiki.org/StockholmTools>

=head1 OUTPUT FORMAT

The output is Stockholm format, interspersed with nonprintable control characters that will be interpreted as colors by an ANSI-compliant terminal (including e.g. VT100 or X11 terminals).
Alternatively, the colors can be encoded in HTML using the C<-html> option.

In addition to the colorizing, some extra rows (summarizing the base-pairing properties) are added to the output alignment.
Each row is prefixed by Stockholm-format tags:

=over 6

=item C<#=GC CS ...>

The characters in this line indicate the number of compensatory substitutions in each column, by counting the number of Watson-Crick or wobble basepairs to be found if that column is considered in tandem with its base-paired partner.

Valid characters include 0, 1, 2, 3, 4, 5 and 6.

=item C<#=GC BP ...>

The characters in this line count the number of rows which have a Watson-Crick or canonical nucleotide pair, in the site of this column and its base-paired partner.

Valid characters include 0-9 (meaning zero to nine rows), a-z (meaning 10 to 35 rows) or + (meaning 36 or more rows).

=item C<#=GC SN ...>

The characters in this line indicate the stem number of each base-paired column.
(Each stem also receives a separate color.)

Valid characters include 0-9 (stems zero through nine), A-Z (stems 10 through 35) or a-z (stems 36 through 61).
In the event of an alignment with more than 62 stems, the numbers wrap round again, so that 0 is also used to denote stems 62, 124, 186 (etc.)

=item C<#=GF CS ...>

This row includes a summary of information about compensatory mutations that may be useful as a quick screen of how much the alignment looks like noncoding RNA.
(See DESCRIPTION for a definition of "compensatory mutations".)

The information includes (i) the total number of columns; (ii) the number of columns that are base-paired; (iii) the number of columns that display compensatory mutations
(see C<-cov> option); and (iv) the number of unbroken stems.

=back

=head1 OPTIONS

=over 12

=item C<-h>

Prints a short help message.

=item C<-less>

Pages the output through C<less -r> (the C<-r> option preserves terminal control characters).

By default, the output is piped to C<cat>.

=item C<-cols N>

Sets the number of columns that will be printed per row of the output.

=item C<-gc TAG>

Allows the user to specify another '#=GC' tag instead of 'SS_CONS' to be used for the consensus secondary structure.

=item C<-cov>

Only colors those columns with compensatory mutations.
(See DESCRIPTION for a definition of "compensatory mutations".)

=item C<-can>

Only colors canonical (Watson-Crick) or wobble base-pairs.

=item C<-inv SEQNAME>

Highlights the row named SEQNAME, by inverting the foreground and background colors.

=item C<-diff SEQNAME>

Highlights positions that are mutated relative to the row named SEQNAME, by inverting the foreground and background colors.

=item C<-ref SEQNAME>

Uses an alternate color scheme, more similar to the sscolor.pl program.
The coloring depends on a reference sequence (SEQNAME).

Sites are colored by whether they are (i) single-stranded or base-paired; (ii) identical or mutated relative to the reference sequence; and (iii) [in the case of base-pairs only] canonical/wobble vs non-canonical.

Specifically, the alternate color scheme is as follows:

=over 2

=item White

Consensus structure can be single-stranded or base-paired.

No substitution relative to reference sequence.

=item Green

Consensus structure implies Watson-Crick or wobble base-pair.

Double substitution relative to reference sequence.

=item Blue

Consensus structure implies Watson-Crick or wobble base-pair.

Single substitution relative to reference sequence.

=item Red

Consensus structure implies non-canonical base-pair.

=item Magenta

Consensus structure implies position is single-stranded.

Single substitution relative to reference sequence.

=back

The C<-ref> option cannot be used with the options for the primary color scheme
(C<-cov>, C<-can>, C<-diff> or C<-inv>).

More info on the sscolor.pl tool can be found here:

L<http://biowiki.org/RnaAlignmentViewers>

=item C<-summ>

Instead of printing the alignment, this just displays the information about compensatory mutations that is normally included in the C<#=GF CS> line.

The information includes (i) the total number of columns; (ii) the number of columns that are base-paired; (iii) the number of columns that display compensatory mutations
(see C<-cov> option); and (iv) the number of unbroken stems.

See DESCRIPTION for a definition of "compensatory mutations".

=item C<-html>

Prints the output in HTML instead of using ANSI terminal control codes.

The colors are defined in a style sheet, for easy editing.

=back

=head1 AUTHOR

Ian Holmes - <ihh@berkeley.edu>

=head1 LICENSE

This code is released under the GNU Public License v3.0.

=cut

sub make_header {
    my $header = <<END;
<?xml version="1.0" encoding="iso-8859-1"?>
	<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US"><head><title>$stockfile</title>
<style type="text/css">
<!--/* <![CDATA[ */
.R {
  background-color: #000000;  /* black */
  color: #FF0000;    /* red */
}
.G {
  background-color: #000000;  /* black */
  color: #00FF00;    /* green */
}
.Y {
  background-color: #000000;  /* black */
  color: #FFFF00;    /* yellow */
}
.B {
  background-color: #000000;  /* black */
  color: #0000FF;    /* blue */
}
.M {
  background-color: #000000;  /* black */
  color: #FF00FF;    /* magenta */
}
.C {
  background-color: #000000;  /* black */
  color: #00FFFF;    /* cyan */
}
.W {
  background-color: #000000;  /* black */
  color: #FFFFFF;    /* white */
}

.IR {
  color: #000000;    /* black */
  background-color: #FF0000;  /* red */
}
.IG {
  color: #000000;    /* black */
  background-color: #00FF00;  /* green */
}
.IY {
  color: #000000;    /* black */
  background-color: #FFFF00;  /* yellow */
}
.IB {
  color: #000000;    /* black */
  background-color: #0000FF;  /* blue */
}
.IM {
  color: #000000;    /* black */
  background-color: #FF00FF;  /* magenta */
}
.IC {
  color: #000000;    /* black */
  background-color: #00FFFF;  /* cyan */
}
.IW {
  color: #000000;    /* black */
  background-color: #FFFFFF;  /* white */
}

/* ]]> */-->
</style>
</head><body bgcolor="#000000">
<pre>
END
    $header .= '<span class="W">';
return $header;
}
