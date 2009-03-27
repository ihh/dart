#!/usr/bin/perl -w

use CGI qw/:standard/;
use Digest::MD5 qw(md5_base64);

# Known problems:
# - need to re-think the whole filestore thing, current system is v scrappy
# - looks like carriage returns (chr13) *still* aren't stripped from uploads
# - synchronous behavior is v dumb, not very useable


# external binaries
my $qsub = "qsub -cwd -v PATH -b y -e /dev/null -sync y";
$qsub = undef;  # disables queueing; comment out this line to use qsub
my $DARTDIR = exists($ENV{'DARTDIR'}) ? $ENV{'DARTDIR'} : "/Users/yam/dart";
my $xrateExecutable = defined($DARTDIR) ? "$DARTDIR/bin/xrate" : "xrate";

# webserver paths
$WEBDOCS = "/Library/Webserver/Documents/xratecgi";  # somewhere under DocumentRoot or an Alias in httpd.conf
$WEBROOT = "http://localhost/xratecgi";

# paths to the store
ensureDirExists ($ADIR = "align");
ensureDirExists ($GDIR = "grammar");
ensureDirExists ($CLADOGRAM = "cladogram");
ensureDirExists ($BUBBLEPLOT = "bubbleplot");

$TMP = "tmp";
$SAVE = "saved";

# temp file path
$TMPPATH = "/tmp";

# get CGI params
my $task = param('task');
my $alignment = param('alignment');
my $grammar = param('grammar');
my $async = param('async');
my $alignLoadPath = safeParam('alignLoadPath');
my $alignSavePath = safeParam('alignSavePath');
my $gramLoadPath = safeParam('gramLoadPath');
my $gramSavePath = safeParam('gramSavePath');
my $alignUpload = upload('alignUpload');
my $gramUpload = upload('gramUpload');

$task = "" unless defined $task;
$alignment = "" unless defined $alignment;
$grammar = "" unless defined $grammar;

# make temporary files & MD5 hashes
($TMPALIGN, $MD5ALIGN) = tmpfile ($alignment);
($TMPGRAMMAR, $MD5GRAMMAR) = tmpfile ($grammar);
$MD5 = my_md5 ($alignment . $grammar);

# path suffices derived from CGI params
$TMPMD5 = "$TMP/$MD5";

$ATMP = "$ADIR/$TMPMD5";
$AFETCH = "$ADIR/$alignLoadPath";
$ASTORE = "$ADIR/$SAVE/$alignSavePath";

$GTMP = "$GDIR/$TMPMD5";
$GFETCH = "$GDIR/$gramLoadPath";
$GSTORE = "$GDIR/$SAVE/$gramSavePath";

$TREE = "$CLADOGRAM/$MD5ALIGN";
$BUBBLE = "$BUBBLEPLOT/$MD5GRAMMAR";

# other vars
my $xrate = "$xrateExecutable $TMPALIGN";  # xrate command line
my $watch;  # where the results of the current task can be found


# perform task
if (length $task) {

    if ($task eq "tree") {
	watch ($ATMP);
	qsub ("$xrate -e $TMPGRAMMAR", "$WEBDOCS/$ATMP");
	$alignment = results();


    } elsif ($task eq "markup") {
	watch ($ATMP);
	qsub ("$xrate -g $TMPGRAMMAR", "$WEBDOCS/$ATMP");
	$alignment = results();


    } elsif ($task eq "train") {
	watch ($GTMP);
	qsub ("$xrate -g $TMPGRAMMAR -t $WEBDOCS/$GTMP");
	$grammar = results();
	param ('gramSavePath', $TMPMD5);


    } elsif ($task eq "loadAlignment") {
	watch ($AFETCH);
	$alignment = results();
	param ('alignSavePath',"");

    } elsif ($task eq "saveAlignment") {
	store ($alignment, $ASTORE);

    } elsif ($task eq "uploadAlignment") {
	$alignment = "";
	while (<$alignUpload>) {
	    $alignment .= $_;
	}
	$alignment =~ s/@{[chr(13)]}//g;   # remove carriage returns
	param('alignSavePath',param('alignUpload'));  # set the save name


    } elsif ($task eq "loadGrammar") {
	watch ($GFETCH);
	$grammar = results();
	param ('gramSavePath',"");

    } elsif ($task eq "saveGrammar") {
	store ($grammar, $GSTORE);

    } elsif ($task eq "uploadGrammar") {
	$grammar = "";
	while (<$gramUpload>) {
	    $grammar .= $_;
	}
	$grammar =~ s/@{[chr(13)]}//g;   # remove carriage returns
	param('gramSavePath',param('gramUpload'));  # set the save name


    } elsif ($task eq "drawtree") {
	watch ($TREE);
#     (stop if file $TREE already exists)
#    renderTree $TMPALIGN >$TREE
#     (Phylodendron, drawphyl, drawpstree or equivalent)

    } elsif ($task eq "drawbubbles") {
	watch ($BUBBLE);
#     (stop if file $BUBBLE already exists)
#    renderBubbleplots $TMPGRAMMAR >$BUBBLE
#     (interface to Rob's bubbleplot modules)
    }
}


# if running in "async" mode, notify client that results are ready
if (defined $async) {
    pingClient();

} else {
# not running in "async" mode, so print CGI form

    print
	header,
	start_html('xrate CGI'),
	h1('xrate CGI'),

#	"task=$task",

	start_multipart_form (-name=>'form'),

	hidden(-name=>'task',
	       -default=>'null'),

	hr,
	div,
	em('Alignment'),

	br,
	'<textarea name="alignment" rows="10" cols="100">',
	$alignment,
	'</textarea>',

	br,
	'Save alignment as ',
	textfield(-name=>'alignSavePath',
		  -onFocus=>"task.value='saveAlignment'"),

	button(-name=>'Save',
	       -onClick=>"task.value='saveAlignment';submit()"),

	br,
	'Load alignment ',
	popup_menu(-name=>'alignLoadPath',
		   -values=>[savedAlignments()],
		   -onChange=>"task.value='loadAlignment';submit()"),

	' or ',
	filefield(-name=>'alignUpload',
		  -default=>'',
		  -size=>50,
		  -maxlength=>80,
		  -onChange=>"task.value='uploadAlignment';submit()"),,


	hr,
	div,
	em('Grammar'),

	br,
	'<textarea name="grammar" rows="10" cols="100">',
	$grammar,
	'</textarea>',

	br,
	'Save grammar as ',
	textfield(-name=>'gramSavePath',
		  -onFocus=>"task.value='saveGrammar'"),

	button(-name=>'Save',
	       -onClick=>"task.value='saveGrammar';submit()"),

	br,
	'Load grammar ',
	popup_menu(-name=>'gramLoadPath',
		   -values=>[savedGrammars()],
		   -onChange=>"task.value='loadGrammar';submit()"),

	' or ',
	filefield(-name=>'gramUpload',
		  -default=>'',
		  -size=>50,
		  -maxlength=>80,
		  -onChange=>"task.value='uploadGrammar';submit()"),,



	hr,
	div,

	button(-name=>'Build Tree',
	       -onClick=>"task.value='tree';submit()"),

	button(-name=>'Annotate Alignment',
	       -onClick=>"task.value='markup';submit()"),

	button(-name=>'Train Grammar',
	       -onClick=>"task.value='train';submit()"),

	end_form,
	hr,

	end_html;
}

# shell metacharacter-safe CGI::param wrapper
sub safeParam {
    my $param = shift;
    my $val = param ($param);
    $val = "" unless defined $val;
    $val =~ s/[^A-Za-z0-9_\-\+\.\:\(\)\/ ]//g;  # eliminate unsafe characters
    return $val;
}

# subroutine to dump data into a temporary file & return filename
sub tmpfile {
    my $data = shift;
    return ("", "") unless defined($data) && length($data);
    my $md5 = my_md5 ($data);
    my $file = $TMPPATH.'/'.$md5;
    unless (-e $file) {
	local *FILE;
	open FILE, ">$file" or die "Couldn't create tmpfile $file: $!";
	binmode FILE;
	print FILE $data;
	close FILE;
    }
    return ($file, $md5);
}

# subroutine to run qsub; returns job ID, or undef if queueing disabled
sub qsub {
    my ($command, $outfile) = @_;
    $outfile = "/dev/null" unless defined $outfile;
    if (defined $qsub) {
	system "$qsub -o $outfile $command";
    } else {
	system "$command >$outfile";
    }
}

# subroutine stub for notifying the client where to look asynchronously for the results of the current task
sub watch {
    my ($path) = @_;
    # ensure directory path exists
    ensurePathExists ($path);
    if (defined $async) {
	print header('text/plain'), "$WEBROOT/$path", "\n";
    }
    $watch = $path;
}

# subroutine to ensure all intermediate directories in a path to a file exist
sub ensurePathExists {
    my ($path) = @_;
    my @dir = split /\/+/, $path;
    for (my $i = 0; $i < @dir-1; ++$i) {
	my $subpath = "$WEBDOCS/" . join ('/', @dir[0..$i]);
	last if -e $subpath && !-d $subpath;  # fails silently if can't create an intermediate directory
	mkdir $subpath or die $! unless -e $subpath;
    }
}

# hacky variant of ensurePathExists for main directories
sub ensureDirExists {
    my ($path) = @_;
    ensurePathExists ("$path/.");  # hack
}

# subroutine to get the results of the current task
sub results {
    return `cat $WEBDOCS/$watch`;
}

# subroutine stub for asynchronously notifying the client that the current task is complete
sub pingClient {
    # do something here, e.g. XML message to client containing ($task,$watch)
}

# subroutine to store a file
sub store {
    my ($data, $filename) = @_;
    ensurePathExists ($filename);
    my $path = "$WEBDOCS/$filename";
    local *FILE;
    open FILE, ">$path";
    print FILE $data;
    close FILE;
}

# subroutine to return list of files saved on server
sub listDir {
    my $root = shift;
    my @align = grep defined() && /\S/, `find $root -type f -print`;
    grep chomp, @align;
    grep s/\s*$//, @align;
    grep s!^$root/*!!g, @align;
    @align = grep !/^$TMP/, @align;
    return @align;
}

# subroutine to return list of alignments saved on server
sub savedAlignments {
    return listDir ("$WEBDOCS/$ADIR");
}

# subroutine to return list of grammars saved on server
sub savedGrammars {
    return listDir ("$WEBDOCS/$GDIR");
}

# md5_base64 wrapper
sub my_md5 {
    my ($data) = @_;
    my $md5 = md5_base64 ($data);
    $md5 =~ s!/!!g;
    return $md5;
}
