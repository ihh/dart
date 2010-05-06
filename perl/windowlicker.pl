#!/usr/bin/perl -w

use strict;
use warnings;

use FileHandle;
use IPC::Open3;
use IO::Select;
use Symbol;

use Stockholm;
use Stockholm::Database;

# constants
my $DARTDIR = 'DARTDIR';  # name of DART's environment variable; windowlicker will look for xrate in $ENV{$DARTDIR}/bin
my $windowLickerTag = "WL";  # Stockholm '#=GF' tag for windowlicker info
my $ID = "ID";  # '#=GC ID' line
my $originalID = "ORIGINAL_ID";  # '#=GC ORIGINAL_ID' line

# regular expressions to filter out of the xrate stderr stream
# be careful to escape regexp special characters
# this is all a bit fragile and hacky...
my @xrateErrRegexp = ("waiting for alignments on standard input",
		      "Warning: duplicate sequence name");
my $xrateErrRegexp = "(" . join ("|", @xrateErrRegexp) . ")";  # combined pattern

my $blockSize = 4096;  # max number of bytes to read in any one chunk from the xrate stdout or stderr streams

# in the GFF output, we want to be GFF3-compliant by default, which means using SOFA terms for the "type" field.
# "region" is the most generic SOFA term for an interior, so we use that as the default "type".
# See Eilbeck et al, 2005. "The Sequence Ontology: a tool for the unification of genome annotations." Genome Biology 6:5.
my $metadataType = "region";      # default GFF "type" field for Stockholm metadata
my $windowNamePrefix = "WINDOW:"; # prefix for name of sliding-window sub-alignment

# get program name
my $progname = $0;
$progname =~ s!.*/!!;

# command-line options
# see $usage string (below) for documentation of these
my $start = 1;
my $end = undef;
my $windowSize = 200;
my $delta;

my $refSeqName;
my $refSeqArgs = "";
my $maxGapFraction = .8;
my $minEntropy = 2;

my $strand = 'f';
my $oneShot = 0;

my $terse = 0;
my $cat = 0;
my $verbose = 0;

my $chunkSize = 100;

my ($gff_file, $gffRefSeq, $wig_file);
my $filename;

# guess path to xrate
my $xrate = "";
$xrate = $ENV{$DARTDIR}.'/bin/' if defined $ENV{$DARTDIR};
$xrate .= "xrate";

# usage text
my $usage;
$usage  = "Usage: $progname [arguments] [Stockholm alignment file] [-- <xrate-args>]\n";
$usage .= "      -s <n>  Start column (default is first column in alignment, i.e. 1)\n";
$usage .= "      -e <n>  End column (inclusive; default is last column in alignment)\n";
$usage .= "      -w <n>  Window size (default is $windowSize)\n";
$usage .= "      -d <n>  \"Delta\" = step size (default is half of window size)\n";
$usage .= "  -dir <dir>  strand: 'f' or 'r', or 'fr' for both (default is '$strand')\n";
$usage .= "   -r <name>  reference sequence in input alignment\n";
$usage .= "      -g <n>  maximum gap fraction in reference sequence (default is $maxGapFraction)\n";
$usage .= "      -b <n>  minimum bits per dinucleotide of reference sequence (default is $minEntropy)\n";
$usage .= "          -o  One-shot mode: overrides -w, sets window size = end + 1 - start\n";
$usage .= "   -x <path>  path to xrate executable (default is '$xrate')\n";
$usage .= "      -c <n>  number of windows to run per xrate process (default is $chunkSize)\n";
$usage .= "        -cat  don't pipe through xrate, just cat windowed alignments to stdout\n";
$usage .= "      -terse  don't print sequence data, just Stockholm metadata\n";
$usage .= "          -v  turn on verbose (log more stuff to stderr, default is off)\n";
$usage .= " -gff <file>  copy GFF annotations to separate file, with Stockholm metadata\n";
$usage .= "  -gr <name>  reference sequence in 1st column of GFF file (default is input filename)\n";
$usage .= " -wig <file>  copy Wiggle annotations to separate file\n";
$usage .= "\n";
$usage .= "Use '-' for alignment filename to accept alignments on standard input\n";
$usage .= " (this can be omitted, unless you also want to specify arguments to xrate)\n";
$usage .= "\n";

# parse command-line opts
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-./) {
	if ($arg eq '-s') {
	    defined ($start = shift) or die $usage;
	} elsif ($arg eq '-e') {
	    defined ($end = shift) or die $usage;
	} elsif ($arg eq '-w') {
	    defined ($windowSize = shift) or die $usage;
	} elsif ($arg eq '-o') {
	    $oneShot = 1;
	} elsif ($arg eq '-d') {
	    defined ($delta = shift) or die $usage;
	} elsif ($arg eq '-x') {
	    defined ($xrate = shift) or die $usage;
	} elsif ($arg eq '-c') {
	    defined ($chunkSize = shift) or die $usage;
	} elsif ($arg eq '-r') {
	    defined ($refSeqName = shift) or die $usage;
	} elsif ($arg eq '-g') {
	    defined ($maxGapFraction = shift) or die $usage;
	    $refSeqArgs .= "$arg $maxGapFraction ";
	} elsif ($arg eq '-b') {
	    defined ($minEntropy = shift) or die $usage;
	    $refSeqArgs .= "$arg $minEntropy ";
	} elsif ($arg eq '-dir') {
	    defined ($strand = shift) or die $usage;
	} elsif ($arg eq '-terse') {
	    $terse = 1;
	} elsif ($arg eq '-gff') {
	    defined ($gff_file = shift) or die $usage;
	} elsif ($arg eq '-gr') {
	    defined ($gffRefSeq = shift) or die $usage;
	} elsif ($arg eq '-wig') {
	    defined ($wig_file = shift) or die $usage;
	} elsif ($arg eq '-cat') {
	    $cat = 1;
	} elsif ($arg eq '-v') {
	    $verbose = 1;
	} elsif ($arg eq '--') {
	    last;
	} else {
	    die $usage, "Unknown option: $arg\n";
	}
    } else {
	die $usage if defined $filename;
	$filename = $arg;
    }
}

if (length($refSeqArgs) && !defined($refSeqName)) {
    die $usage, "Reference sequence needed for arguments: $refSeqArgs\n";
}

if (!defined $filename) {
    warn "[waiting for alignments on standard input]\n";
    $filename = '-';
}

# get xrate args
my @xrateArgs = @ARGV;

# read alignment
my $stockdb = Stockholm::Database->from_file ($filename, $verbose);
die "No alignments in input\n" if @$stockdb < 1;
die "Can't process more than one input alignment\n" if @$stockdb > 1;
my $stock = $stockdb->[0];
die "Empty alignment\n" unless $stock->columns;

# post-process command line args, in light of alignment
$end = $stock->columns unless defined $end;
die "Nonsense argument: end ($end) > alignment length (", $stock->columns, ")\n" if $end > $stock->columns;

# window statistics
$windowSize = $end + 1 - $start if $oneShot;

$delta = int ($windowSize / 2) unless defined $delta;
die "Nonsense arguments: delta ($delta) > window size ($windowSize)\n" if $delta > $windowSize;

my $overlap = $windowSize - $delta;

# reference sequence
my $refSeq;
if (defined $refSeqName) {
    die "Reference sequence '$refSeqName' not found\n" unless exists $stock->seqdata->{$refSeqName};
    $refSeq = $stock->seqdata->{$refSeqName};
}

# GFF (and WIG) reference sequence
unless (defined $gffRefSeq) {
    if (defined $refSeqName) {
	# use reference sequence name, if defined
	$gffRefSeq = $refSeqName;

    } elsif (defined $stock->get_gf ($ID)) {
	# use original '#=GF ID' as ref seq name, if there was one
	$gffRefSeq = $stock->get_gf ($ID);

    } else {
	# use input file name as ref seq name
	# (if input is stdin, replace the '-' because it looks
	# confusing in a GFF file)
	$gffRefSeq = ($filename eq '-') ? 'stdin' : $filename;
    }
}

# translate window strand from command-line syntax ("f" or "r", etc.) into "+" or "-" for Stockholm.pm
$strand = lc $strand;
$strand =~ tr/fp5bmnr3/+++\-\-\-\-\-/;
die "Don't understand '$strand' as a strand direction\n" unless grep ($strand eq $_, qw(+ - +- -+));

# autoflush on
$| = 1;

# open GFF file for annotations and #=GF,#=GC metadata
local *GFF_FILE;
if (defined $gff_file) {
    open GFF_FILE, ">$gff_file" or die "Couldn't open '$gff_file': $!\n";
    select((select(GFF_FILE), $| = 1)[0]);  # autoflush
    print GFF_FILE "##gff-version   3\n";
}

# create structure to hold Wiggle annotations
my %wiggleTrack;

# globals that get used directly by subs
my ($windowsAttempted, $windowsDone) = (0, 0);
my $stockDbIn = Stockholm::Database->new;  # aggregates chunks of windows for xrate
my @stockDbStrands;  # have to keep track of f/r strands for each window in the chunk,
                     # can't use Stockholm::Database for it

# main window loop
my $inLastWin = 0;
window:
for (my $winStart = $start; !$inLastWin; $winStart += $delta) {

    # get window coords
    my $winEnd = $winStart + $windowSize - 1;
    if ($winEnd >= $end) {
	$inLastWin = 1;
	$winEnd = $end;
    }
    my $winLen = $winEnd + 1 - $winStart;

    # skip if refseq is gappy or low-complexity
    if (defined $refSeqName) {
	my $refSubSeq = substr ($refSeq, $winStart - 1, $winLen);
	my $refGaps = $refSubSeq;
	$refSubSeq =~ s/[\-\._]//g;
	$refGaps =~ s/[^\-\._]//g;

	my $gapFrac = length($refGaps) / $winLen;
	if ($gapFrac > $maxGapFraction) {
	    my $gapPercentage = int (100 * $gapFrac + .5);
	    warn "[skipping window from columns $winStart to $winEnd: reference sequence '$refSeqName' is $gapPercentage\% gappy]\n";
	    if (!$inLastWin) {
		next window;
	    } else {
		# no more windows left; need to run chunk
		run_chunk() if @stockDbStrands > 0;
		last window;
	    }
	}

	my ($n, %dinuc);
	for ($n = 0; $n < length($refSubSeq) - 1; ++$n) {
	    ++$dinuc{lc (substr($refSubSeq, $n, 2))};
	}
	my $entropy = 0;
	while (my ($dinuc, $count) = each %dinuc) {
	    my $p = $count / $n;
	    $entropy -= $p * log($p);
	}
	$entropy /= log(2);
	if ($entropy < $minEntropy) {
	    warn "[skipping window from columns $winStart to $winEnd: reference sequence '$refSeqName' has only $entropy bits per dinucleotide]\n";
	    if (!$inLastWin) {
		next window;
	    } else {
		# no more windows left; need to run chunk
		run_chunk() if @stockDbStrands > 0;
		last window;
	    }
	}
    }

    # loop over strands
    for my $winStrand (split (//, $strand)) {

	# translate window strand back into command-line syntax ("f" or "r")
	my $winStrandChar = $winStrand;
	$winStrandChar =~ tr/+-/fr/;

	# get the subalignment
	my $subalign = $stock->subalign ($winStart - 1, $winLen, $winStrand);

	# set the name
	my $windowName = $windowNamePrefix . $winStart . "-" . $winEnd;

	# save window name/identifier as a '#=GF ID' annotation,
	# moving the original annotation (if there was one) to '#=GF ORIGINAL_ID'
	if (defined $subalign->get_gf ($ID)) {
	    $subalign->add_gf ($originalID, $subalign->get_gf ($ID));
	    $subalign->set_gf ($ID, $windowName);
	}
	else {
	    $subalign->add_gf ($ID, $windowName);
	}

	# add windowlicker tag
	$subalign->add_gf ($windowLickerTag, "$progname -o -s $winStart -e $winEnd -dir $winStrandChar");

	# print to stdout or pipe to xrate?
	if ($cat) {
	    print $subalign->to_string ('NOSEQDATA' => $terse);

	} else {
	    # add this window to current chunk of windows for xrate
	    $stockDbIn->add_alignment ($subalign, $verbose);
	    push (@stockDbStrands, $winStrand);
	    $windowsAttempted++;

	    if ((@$stockDbIn == $chunkSize) or $inLastWin) {
		# send chunk of windows to xrate
		run_chunk();
	    }
	}
    }
}

# close GFF file
if (defined $gff_file) {
    close GFF_FILE or die "Couldn't close '$gff_file': $!\n";
}

# write Wiggle file
if (defined $wig_file) {
    local *WIG_FILE;
    open WIG_FILE, ">$wig_file" or die "Couldn't open '$wig_file': $!\n";
    select((select(WIG_FILE), $| = 1)[0]);  # autoflush
    while (my ($track, $data) = each %wiggleTrack) {
	print WIG_FILE map ("$_\n", $track, "fixedStep chrom=$gffRefSeq start=0 step=1", map (ref($_) ? mean($_) : $_, @$data));
    }
    close WIG_FILE or die "Couldn't close '$wig_file': $!\n";
}

# finish
unless ($cat) {
    warn "attempted $windowsAttempted windows on xrate, successfully ran $windowsDone windows\n";
}

exit 0;


################################################################################
#  helper subs
################################################################################

# escape (by converting to hex, a la URL escaping conventions)
# characters that are not allowed in the GFF "attributes" column
sub escapeAttr {
    my ($string) = @_ or die "Invalid use of escapeAttr(string) sub!";
    $string =~ s/([,=;])/uc sprintf("%%%02x",ord($1))/eg;
    return $string;
}

# mean of a list
sub mean {
    my ($listRef) = @_;
    return undef unless @$listRef > 0;
    my $total = 0;
    for my $x (@$listRef) { $total += $x }
    return $total / @$listRef;
}

# run a set of Stockholm alignments (windows) through an instance of xrate
sub run_xrate {
    my $command = "$xrate @xrateArgs";

    # open a tridirectional pipe thingy
    my ($XrateReader, $XrateWriter, $XrateError);
    $XrateError = gensym();
    my $pid = open3 ($XrateWriter, $XrateReader, $XrateError, $command);
    $SIG{'PIPE'} = sub { die "xrate failed (check command-line options? they were: $command)\n" };

# leave for debugging (might be suffering from buffering...)
#    warn "xrate:\nstrace -p $pid\nwindowlicker:\nstrace -p $$\n";
#    sleep 10;

    # send the alignments to xrate
    print $XrateWriter $stockDbIn->to_string;
    close $XrateWriter;

    # NB: the following code (and comments) for reading from pipes
    #     is thanks to http://www.perlmonks.org/?node_id=151886

    my $sel = new IO::Select ($XrateReader, $XrateError);
    my ($xrateOut, $xrateErr) = ('', '');

    # $sel->can_read will block until there is data available on one or more fhs
    while (my @ready = $sel->can_read) {

	# now we have a list of all fhs that we can read from
	foreach my $fh (@ready) { # loop through them

	    my $line;
	    # read up to $blockSize bytes from this fh.
	    # if there is less than $blockSize bytes, we'll only get
	    # those available bytes and won't block.  If there
	    # is more than $blockSize bytes, we'll only read $blockSize and
	    # wait for the next iteration through the loop to
	    # read the rest.
	    my $len = sysread $fh, $line, $blockSize;
	    if (not defined $len) {
		# There was an error reading
		die "Error from child: $!\n";

	    } elsif ($len == 0){
		# Finished reading from this FH because we read
		# 0 bytes.  Remove this handle from $sel.
		# we will exit the loop once we remove all file
		# handles ($outfh and $errfh).
		$sel->remove ($fh);
		next;

	    } else { # we read data alright
		#print "Read $len bytes from $fh\n";
		if ($fh == $XrateReader) {
		    $xrateOut .= $line;

		} elsif ($fh == $XrateError) {
		    $xrateErr .= $line;

		} else {
		    die "Shouldn't be here";
		}
	    }
	}
    }

    # read xrate's stdout into a Stockholm database
    my $stockDbOut = Stockholm::Database->from_string ($xrateOut, $verbose);

    # filter out xrate's "[waiting for alignments on stdin]" and "duplicate sequence name" warnings
    # (HACKY! could use a temp file instead...)
    my @xrate_err = split ("\n", $xrateErr);
    @xrate_err = map {/$xrateErrRegexp/ ? () : "[stderr] $_\n"} @xrate_err;   # remove standard crap from xrate stderr
    warn @xrate_err if @xrate_err;

    # reap child xrate process, so there are no zombies hanging around
    waitpid ($pid, 0);

    # if xrate failed, we probably don't want to waste CPU time running the rest of
    # the segment until the failure is fixed (odds are the failure will just happen
    # again on a later window chunk), so quit with error
    die "xrate returned with nonzero exit status (invoked using: $command)\n" unless $? == 0;

    return $stockDbOut;
}

# process and output a chunk of windows after they are run through xrate
sub run_chunk {
    die 'this should not occur' unless (@stockDbStrands and (@stockDbStrands == @$stockDbIn));

    my $stockDbOut = run_xrate();

    if (@$stockDbIn != @$stockDbOut) {
	warn "WARNING: xrate lost ", scalar (@$stockDbIn) - scalar (@$stockDbOut), " windows (out of ", scalar (@$stockDbIn), ") when processing chunk.\n";
    }

    # go through windows in this chunk, parse and output them
    for (my $i = 0; $i < @$stockDbOut; $i++) {
	my $winIn = $stockDbIn->[$i];
	my $winOut = $stockDbOut->[$i];
	my $winStrand = $stockDbStrands[$i];
	my $windowName = $winIn->get_gf ($ID);
	$windowName =~ /^$windowNamePrefix(\d+)-(\d+)$/ or die 'This should not occur!';
	my ($winStart, $winEnd) = ($1, $2);

	# if xrate output parsed OK, add coords to alignment & print
	if (defined $winOut->columns) {
	    if (defined $gff_file) {
		# attributes for window metadata GFF line
		my @tag_value = ("ID=" . escapeAttr ($windowName) . ";Note=window metadata");

		# add new #=GF tags
		my @gff;
		for my $gf_tag (sort keys %{$winOut->gf}) {
		    if ($gf_tag eq "GFF") {
			# handle GFF tags separately; correct the start & end coords and point seqname to the main alignment
			my %old_val = map (($_ => 1), @{$winIn->gf_($gf_tag)});
			my @new_val = @{$winOut->gf_($gf_tag)};
			for my $gff (@new_val) {
			    unless (exists $old_val{$gff}) {
				my @gff_field = split (/\t/, $gff, 9);
				$gff_field[0] = $gffRefSeq;
				$gff_field[3] += $winStart - 1;
				$gff_field[4] += $winStart - 1;
				my $newAttributes = "Parent=" . escapeAttr ($windowName);
				$gff_field[8] = ($gff_field[8] =~ /^\s*\.?\s*$/ ?
						 $newAttributes :
						 "$gff_field[8];$newAttributes");
				push @gff, join ("\t", @gff_field) . "\n";
			    }
			}

		    } elsif ($gf_tag eq "WIG") {
			# handle WIG tags separately: correct the start points
			# makes some strong
			my ($currentWigTrack, $pos);
			for my $wig (@{$winOut->gf_($gf_tag)}) {
			    if ($wig =~ /^track\b/) {
				$wig =~ s/\b(name=\S+)\b/$1$winStrand/;
				$wiggleTrack{$wig} = [] unless exists $wiggleTrack{$wig};
				$currentWigTrack = $wiggleTrack{$wig};
				$pos = $winStrand eq '+' ? ($winStart-1) : ($winEnd-1);
			    } elsif (defined $currentWigTrack) {
				if ($wig =~ /^fixedStep\b/) {
				    # ignore fixedStep lines
				} else {
				    if (defined $currentWigTrack->[$pos]) {
					# store overlapping values for later averaging
					unless (ref $currentWigTrack->[$pos]) {
					    $currentWigTrack->[$pos] = [$currentWigTrack->[$pos]];
					}
					push @{$currentWigTrack->[$pos]}, $wig + 0;
				    } else {
					$currentWigTrack->[$pos] = $wig + 0;
				    }
				    $pos += $winStrand eq '+' ? 1 : -1;
				}
			    }
			}

		    } else {
			my $old_val = join "", @{$winIn->gf_($gf_tag)};
			my $new_val = join "", @{$winOut->gf_($gf_tag)};
			if ($new_val ne $old_val) {
			    $gf_tag = escapeAttr ($gf_tag);
			    $new_val = escapeAttr ($new_val);
			    push @tag_value, "stock_GF_$gf_tag=$new_val";
			}
		    }
		}

		# add new #=GC tags
		for my $gc_tag (sort keys %{$winOut->gc}) {
		    my $old_val = $winIn->gc_($gc_tag);
		    my $new_val = $winOut->gc_($gc_tag);
		    if ($new_val ne $old_val) {
			$gc_tag = escapeAttr ($gc_tag);
			$new_val = escapeAttr ($new_val);
			push @tag_value, "stock_GC_$gc_tag=$new_val";
		    }
		}

		print GFF_FILE @gff;
		print GFF_FILE join ("\t", $gffRefSeq, $progname, $metadataType, $winStart, $winEnd, '.', $winStrand, '.', join (";", @tag_value)), "\n";
	    }

	    print $winOut->to_string ('NOSEQDATA' => $terse);
	    $windowsDone++;

	} else {
	    warn "WARNING: no columns in xrate output alignment, skipping window.\n";
	}
    }  # end loop over windows in the chunk

    $stockDbIn = Stockholm::Database->new;  # reset for next chunk
    @stockDbStrands = ();

    warn "ran chunk OK, it contained ", scalar @$stockDbOut, " windows\n";
}
