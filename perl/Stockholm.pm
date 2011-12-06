#!/usr/bin/perl -w

=head1 NAME

Stockholm.pm

=head1 SYNOPSIS

Lightweight Perl module encapsulating a Stockholm multiple alignment.

For more detail on Stockholm format itself, see the following URL:

http://biowiki.org/StockholmFormat

=head1 METHODS

=cut

package Stockholm;

use strict;
use vars '@ISA';

use Carp qw(carp cluck croak);

# define the gap alphabet
my $gapChars = '-._';
my $gapCharsRegexp = '\-\._';
sub gapCharsRegexp { return $gapCharsRegexp }

=head2 new

    my $stock = Stockholm->new();

Creates an empty Stockholm object.

=cut

# constructor
sub new {
    my ($class) = @_;
    my $self = {
	'seqname' => [],  # order of sequences
	'seqdata' => {},
	'gf' => {},   # Generic File annotation (freeform)
	'gc' => {},   # Generic Consensus annotation (by-column)
	'gs' => {},   # Generic Sequence annotation (by-sequence, freeform)
	'gr' => {},   # Generic by-Row annotation (by-sequence, by-column)
	'gfOrder' => []  # order of GF lines
    };
    bless $self, ref($class) || $class;
    return $self;
}

=head2 from_file

    my $stock = Stockholm->from_file ($filename);

Creates a new Stockholm object and reads it from a file.

=cut

# from_file method
sub from_file {
    my ($class, $filename) = @_;

    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    my $self = $class->from_filehandle (\*FILE);
    close FILE;

    return $self;
}

=head2 from_filehandle

    my $stock = Stockholm->from_filehandle ($filehandle);

Creates a new Stockholm object and reads it from a filehandle.

=cut

# from_filehandle method
sub from_filehandle {
    my ($class, $filehandle) = @_;

    my $self = $class->new;
    local $_;
    while (<$filehandle>) {
	last if $self->parse_input_line ($_);
    }
    carp "Warning: alignment is not flush" unless $self->is_flush();

    return $self;
}

=head2 to_file

    $stock->to_file ($filename);

Writes a Stockholm object to a file.

=cut

# to_file method
sub to_file {
    my ($self, $filename, $maxcols) = @_;

    local *FILE;
    open FILE, ">$filename" or croak "Couldn't open '$filename' for writing: $!";
    print FILE $self->to_string ($maxcols);
    close FILE or croak "Couldn't close '$filename': $!";;

    return $filename;
}

=head2 to_string

    print $stock->to_string ($maxcols)
    print $stock->to_string ($maxcols, ARG1=>VAL1, ARG2=>VAL2 ...)
    print $stock->to_string (MAXCOLS=>$maxcols, ARG1=>VAL1, ARG2=>VAL2 ...)

Returns the object as a Stockholm-formatted string.

ARGs can include...
        MAXCOLS    -- set this to an integer $maxcols to limit the maximum number of displayed screen columns to $maxcols. Alternatively, set this to the string "SCREEN" to use all available screen columns.
        NOSEQDATA  -- set this to 1 to suppress printing the actual sequence data in the alignment (for example if you are only interested in the annotation lines)
        COLOR      -- set this to a subroutine reference that, when called with (residue,rowname,colnum), returns 2-elt list of (fg,bg) colors. Alternatively, set this to the string "AMINO" for amino acid coloring.

=cut

sub to_string {
    my ($self, @args) = @_;
    my (%args, $maxcols);
    if (@args % 2 == 1) {
	$maxcols = shift @args;
	%args = @args;
    } else {
	%args = @args;
	$maxcols = $args{'MAXCOLS'};
    }
    $maxcols = (`tput cols` + 0) || 80 if defined($maxcols) && $maxcols eq "SCREEN";
    $maxcols = 80 unless defined $maxcols;     # default 80-column output

    my $color_sub = $args{'COLOR'};
    if (defined($color_sub) && uc($color_sub) eq 'AMINO') { $color_sub = \&amino_color_scheme }

    # init output array
    my @out;
    push @out, "# STOCKHOLM 1.0";

    # determine alignment columns, legend columns & effective columns per line
    my $acols = $self->columns;
    my $lcols = $self->lcols;
    my $colstep = $maxcols < 1 ? $acols : $maxcols - $lcols - 1;
    $colstep = $maxcols if $colstep < 1;  # protect against negative and 0 colstep...

    # GF lines
    # check for gfOrder (insane, fragile Stockholm line ordering strikes again)
    if (@{$self->gfOrder} == map { (@$_) } values %{$self->gf}) {  # gfOrder same number of lines as #=GF block?
	my %gfCursor = map (($_ => 0), keys %{$self->gf});
	foreach my $feature (@{$self->gfOrder}) {
	    push @out, $self->prettify ($lcols, "#=GF $feature", $self->gf_($feature)->[$gfCursor{$feature}++]);
	}
    } else {
	@{$self->gfOrder} = ();  # gfOrder is useless, so flush it
	foreach my $feature (sort { $a cmp $b } keys %{$self->gf}) {
	    push @out, $self->prettify ($lcols, "#=GF $feature", @{$self->gf_($feature)});
	}
    }

    # GS lines
    my @gs_seqname = @{$self->seqname};
    my %gs_seqname_hash = map (($_=>1), grep (!exists $self->seqdata->{$_}, map (keys(%{$self->gs_($_)}), keys %{$self->gs})));
    push @gs_seqname, keys %gs_seqname_hash;
    foreach my $feature (sort { $a cmp $b } keys %{$self->gs}) {
	my $hash = $self->gs_($feature);
	foreach my $seqname (grep (exists($hash->{$_}), @gs_seqname)) {
	    push @out, $self->prettify ($lcols, "#=GS $seqname $feature", @{$hash->{$seqname}});
	}
    }

    my @gcfeat = sort { $a cmp $b } keys %{$self->gc};
    my @gr_seqname = @{$self->seqname};
    my %gr_seqname_hash = map (($_=>1), grep (!exists $self->seqdata->{$_}, map (keys(%{$self->gr_($_)}), keys %{$self->gr})));
    push @gr_seqname, keys %gr_seqname_hash;

    # Loop over columns
    for (my $col = 0; $col < $acols; $col += $colstep) {

	# GC lines
	foreach my $feature (@gcfeat) {
	    push @out, $self->prettify ($lcols, "#=GC $feature",
					safe_substr ($self->gc_($feature), $col, $colstep));
	}
	for (my $i = 0; $i < @gr_seqname; ++$i) {
	    my $seqname = $gr_seqname[$i];
	    # Sequences
#	    warn "Writing cols $col+$colstep for $seqname";
	    if (exists $self->seqdata->{$seqname} && !$args{'NOSEQDATA'}) {
		my $substr = safe_substr ($self->seqdata->{$seqname}, $col, $colstep);
		if (defined $color_sub) {
		    my @substr = split //, $substr;
		    $substr = join ("", map (color_code(&$color_sub($substr[$_],$seqname,$col+$_)).$substr[$_], 0..@substr-1), white_color_code());
		}
		push @out, $self->prettify ($lcols, $seqname, $substr);
	    }
	    # GR lines
	    foreach my $feature (grep (exists ($self->gr->{$_}->{$seqname}), keys %{$self->gr})) {
#		warn "Writing cols $col+$colstep for $seqname $feature";
		push @out, $self->prettify ($lcols, "#=GR $seqname $feature",
					    safe_substr ($self->gr_($feature)->{$seqname}, $col, $colstep));
	    }
	}
	push @out, "";
    }

    # alignment separator
    push @out, "//";

    # convert output array to string & return
    return join ("", map ("$_\n", @out));
}

# ANSI control characters
sub esc { join ("", map (chr(27)."[$_"."m", @_)) }
sub color_code {
    my ($bg, $fg) = @_;
    return esc(27,30+$fg,40+$bg);
}
sub white_color_code { color_code(0,7) }

# color schemes
my %aa_color = ('a' => [7,0],
		'r' => [4,6],
		'n' => [6,7],
		'd' => [1,7],
		'c' => [3,7],
		'q' => [6,0],
		'e' => [1,0],
		'g' => [7,2],
		'h' => [4,7],
		'i' => [2,0],
		'k' => [2,3],
		'l' => [4,0],
		'm' => [3,0],
		'f' => [4,2],
		'p' => [5,0],
		's' => [3,1],
		't' => [3,2],
		'w' => [5,7],
		'y' => [4,0],
		'v' => [2,7]);

sub amino_color_scheme {
    my ($residue, $seqname, $pos) = @_;
    my $cols = $aa_color{lc $residue};
    return defined($cols) ? @$cols : (0,7);
}

=head2 copy

    my $newStock = $stock->copy();

Does a deep-copy, duplicating all information.

=cut

# deep-copy constructor
sub copy {
    my ($self) = @_;
    my $stock = Stockholm->new;

    # Sequence names & data
    @{$stock->seqname} = @{$self->seqname};
    %{$stock->seqdata} = %{$self->seqdata};

    #=GF
    while (my ($feature, $arrayRef) = each %{$self->gf}) {
	$stock->gf->{$feature} = [@$arrayRef];
    }
    @{$stock->gfOrder} = @{$self->gfOrder};

    #=GC
    while (my ($feature, $string) = each %{$self->gc}) {
	$stock->gc->{$feature} = $string;
    }

    #=GR
    while (my ($feature, $seqHash) = each %{$self->gr}) {
	$stock->gr->{$feature} = {%$seqHash};
    }

    #=GS
    while (my ($feature, $seqHash) = each %{$self->gs}) {
	$stock->gs->{$feature} = {%$seqHash};
    }

    # Return
    return $stock;
}

=head2 columns

    my $cols = $stock->columns()

Returns the number of columns in this alignment.

=cut

# Number of columns
sub columns {
    my ($self) = @_;

    return max(map(length($_),
		   values(%{$self->seqdata}),
		   values(%{$self->gc}),
		   map(values(%$_), values(%{$self->gr}))));
}

=head2 sequences

    my $rows = $stock->sequences()

Returns the number of sequences (i.e. rows) in this alignment.

=cut

# Number of sequences (i.e. rows)
sub sequences {
    my ($self) = @_;
    return @{$self->seqname} + 0;
}

=head2 seqname

    my $rowName = $stock->seqname->[$rowIndex];

Returns a reference to an array of sequence names.

=head2 seqdata

    my $row = $stock->seqdata->{$rowName};

Returns a reference to a hash of alignment rows, keyed by sequence name.

=head2 gf

    my @gf = @{$stock->gf->{FEATURE}};
    my @gf = @{$stock->gf_FEATURE};

Returns a reference to an array of all the lines beginning '#=GF FEATURE ...'

=head2 gc

    my $gc = $stock->gc->{FEATURE};
    my $gc = $stock->gc_FEATURE;

Returns the line beginning '#=GC FEATURE ...'

=head2 gs

    my @gs = @{$stock->gs->{FEATURE}->{SEQNAME}};
    my @gs = @{$stock->gs_FEATURE->{SEQNAME}};

Returns a reference to an array of all the lines beginning '#=GS SEQNAME FEATURE ...'

=head2 gr

    my $gr = $stock->gr->{FEATURE}->{SEQNAME};
    my $gr = $stock->gr_FEATURE->{SEQNAME};

Returns the line beginning '#=GR SEQNAME FEATURE ...'

=cut


# catch all methods by default
sub AUTOLOAD {
    my ($self, @args) = @_;
    my $sub = our $AUTOLOAD;
    $sub =~ s/.*:://;

    # check for DESTROY
    return if $sub eq "DESTROY";

    # check for GF, GC, GS, GR tag_feature accessors
    # e.g. $self->GF_ID
    # can also use $self->GF_('ID')
    if ($sub =~ /^(gf|gc|gs|gr)_(\S*)$/i) {
	my ($tag, $feature) = ($1, $2);
	$tag = lc $tag;
	$feature = shift(@args) unless length $feature;
	my $hash = $self->{$tag};
	if (@args) {
	    $hash->{$feature} = shift(@args);   # TODO: check ref-type of arg
	}
	cluck "Warning: ignoring extra arguments to ${tag}_$feature"
	    if @args > 0;
	if (!defined $hash->{$feature}) {
	    $hash->{$feature} = 
		$tag eq "gf" ? [] :
		$tag eq "gc" ? "" :
		$tag eq "gs" ? {} :
		$tag eq "gr" ? {} :
		croak "Unreachable";
	}
	return $hash->{$feature};
    }

    # check for ordinary accessors
    if (exists $self->{$sub}) {
	croak "Usage: $sub() or $sub(newValue)" if @args > 1;
	return
	    @args
	    ? $self->{$sub} = $args[0]
	    : $self->{$sub};
    }

    # croak
    croak "Unsupported method: $sub";
}


# pretty print line(s)
sub prettify {
    my ($self, $lcols, $legend, @data) = @_;
    # This horribly inefficient/redundant series of transformations comes out with something I like (IH, 7/24/07)
    # Trim it down? pah! Like I have nothing better to do
    $legend = sprintf ("% ${lcols}s", $legend);
    $legend =~ s/^(\s+)(\#=\S\S)(\s+\S+)$/$2$1$3/;
    $legend =~ s/^(\s+)(\#=\S\S\s+\S+)(\s+\S+)$/$2$1$3/;
    $legend =~ s/^(\s\s\s\s\s)(\s+)([^\#]\S+)/$1$3$2/;
    return map ("$legend $_", @data);
    # (pukes into cold coffee cup)
}

# legend width (subtract this from maxcols-1 to get number of columns available for sequence display)
sub lcols {
    my ($self) = @_;
    my $lcols = max ($self->maxNameLen,
		     map(length("#=GF $_"), keys(%{$self->gf})),
		     map(length("#=GC $_"), keys(%{$self->gc})));
    while (my ($gr_key, $gr_hash) = each %{$self->gr}) {
	$lcols = max ($lcols, map(length("#=GR $gr_key $_"), keys(%$gr_hash)));
    }
    while (my ($gs_key, $gs_hash) = each %{$self->gs}) {
	$lcols = max ($lcols, map(length("#=GS $gs_key $_"), keys(%$gs_hash)));
    }
    return $lcols;
}

# safer version of substr (doesn't complain if startpoint outside of string)
sub safe_substr {
    my ($string, $start, $len) = @_;
    return "" if !defined($string) || $start > length $string;
    if (defined $len) {
      return substr ($string, $start, $len);
    } else {
      return substr ($string, $start);  # print from $start to end if no length specified
    }
}

=head2 add_gf

    $stock->add_gf (FEATURE, $line)
    $stock->add_gf (FEATURE, $line1, $line2, ...)
    $stock->add_gf (FEATURE, @lines)
    $stock->add_gf (FEATURE, "$line1\n$line2\n$...")

Add '#=GF FEATURE' annotation, preserving the order of the '#=GF' lines
(damn those crazy context-sensitive Stockholm semantics!)

Each list entry gets printed on its own line.

=cut

# Add #=GF annotation, preserving the order of the #=GF lines
# (damn those crazy context-sensitive Stockholm semantics!)
#
# Usage:
#        $stock->add_gf ($gfFeature, $line);
#        $stock->add_gf ($gfFeature, $line1, $line2, ...);
#        $stock->add_gf ($gfFeature, @lines);
#        $stock->add_gf ($gfFeature, "$line1\n$line2\n$...");
#
# Each list entry gets printed on its own line.
#
sub add_gf {
    my ($self, $feature, @data) = @_;
    unless (@data > 0) {
	carp "Missing parameters to 'add_gf'; nothing will be done";
	return;
    }

    foreach my $line (map {$_ eq '' ? '' : split ("\n", $_, -1)} @data) {
	push @{$self->gf_($feature)}, $line;
	push @{$self->gfOrder}, $feature;
    }
}

=head2 set_gf

    $stock->set_gf (FEATURE, $line)
    $stock->set_gf (FEATURE, $line1, $line2, ...)
    $stock->set_gf (FEATURE, @lines)
    $stock->set_gf (FEATURE, "$line1\n$line2\n$...")

Set a '#=GF FEATURE ...' annotation.
The difference between this and 'add_gf' is that here, if this sequence and
feature have an annotation already, it will be overwritten instead of appended to.

As with add_gf, each list entry gets printed on its own line.

=cut

# Set a #=GF annotation.
# The difference between this and 'add_gf' is that here, if this sequence and
# feature have an annotation already, it will be overwritten instead of appended to.
#
# Usage:
#        $stock->set_gf ($gfFeature, $line);
#        $stock->set_gf ($gfFeature, $line1, $line2, ...);
#        $stock->set_gf ($gfFeature, @lines);
#        $stock->add_gf ($gfFeature, "$line1\n$line2\n$...");
#
# Each list entry gets printed on its own line.
sub set_gf {
    my ($self, $feature, @data) = @_;
    unless (defined $self->{gf}->{$feature}) {
	carp "No #=GF feature '$feature' to set, creating new";
	$self->add_gf ($feature, @data);
	return;
    }
    @data = map {$_ eq '' ? '' : split ("\n", $_, -1)} @data;

    ### go to insane lengths to preserve original line ordering

    my %gfCursor = map {$_ => 0} keys %{$self->{gf}};
    my (%new_gfFeature, @new_gfOrder);

    foreach my $curFeat (@{$self->{gfOrder}}) {
	$new_gfFeature{$curFeat} = [] unless defined $new_gfFeature{$curFeat};

	if ($curFeat eq $feature) {
	    if (@data) {
		push (@{$new_gfFeature{$curFeat}}, shift @data);
		push (@new_gfOrder, $curFeat);
	    }
	}
	else {
	    push (@{$new_gfFeature{$curFeat}},
		  $self->{gf}->{$curFeat}->[$gfCursor{$curFeat}++]);
	    push (@new_gfOrder, $curFeat);
	}
    }

    # still have lines left, insert them after last occurence of feature tag
    if (@data) {
	my $insertAt = $#new_gfOrder;
	$insertAt-- while $new_gfOrder[$insertAt] ne $feature;
	splice (@new_gfOrder, ++$insertAt, 0, map {$feature} @data);
	push (@{$new_gfFeature{$feature}}, @data);
    }

    $self->{gf} = {%new_gfFeature};
    $self->{gfOrder} = [@new_gfOrder];
}

=head2 get_gf

    my $gf = $stock->get_gf (FEATURE)

Concatenates and returns all lines beginning '#=GF FEATURE ...'

=cut

# Get a #=GF annotation.
#
sub get_gf {
    my ($self, $feature) = @_;
    if (defined $self->{gf}->{$feature}) {
	return join ("\n", @{$self->{gf}->{$feature}});
    }
    else {
	#carp "No #=GF feature '$feature' in Stockholm object, returning undef";
	return undef;
    }
}

=head2 add_gs

    $stock->add_gs (SEQNAME, FEATURE, $line)
    $stock->add_gs (SEQNAME, FEATURE, $line1, $line2, ...)
    $stock->add_gs (SEQNAME, FEATURE, @lines)
    $stock->add_gs (SEQNAME, FEATURE, "$line1\n$line2\n$...")

Add a GS annotation: '#=GS SEQNAME FEATURE ...'

If such a line already exists for this sequence and feature, append to it.

Each list entry gets printed on its own line.

=cut

# Add a GS annotation: #=GS <sequence> <feature> <free text>
# If one already exists for this sequence and feature, append to it.
#
# Usage:
#        $stock->add_gs ($seq, $gsFeature, $line);
#        $stock->add_gs ($seq, $gsFeature, $line1, $line2, ...);
#        $stock->add_gs ($seq, $gsFeature, @lines);
#        $stock->add_gs ($seq, $gsFeature, "$line1\n$line2\n$...");
#
# Each list entry gets printed on its own line.
#
sub add_gs {
  my ($self, $seq, $feature, @text) = @_;
  unless (@text > 0) {
      carp "Missing parameters to 'add_gs' method; nothing will be done";
      return;
  }
  @text = map {$_ eq '' ? '' : split ("\n", $_, -1)} @text;
  push (@{$self->{gs}->{$feature}->{$seq}}, @text);
}

=head2 set_gs

    $stock->set_gs (SEQNAME, FEATURE, $line)
    $stock->set_gs (SEQNAME, FEATURE, $line1, $line2, ...)
    $stock->set_gs (SEQNAME, FEATURE, @lines)
    $stock->set_gs (SEQNAME, FEATURE, "$line1\n$line2\n$...")

Set a GS annotation: '#=GS SEQNAME FEATURE ...'

The difference between this and 'add_gs' is that here, if this sequence and
feature have an annotation already, it will be overwritten instead of appended to.

Each list entry gets printed on its own line.

=cut

# Set a #=GS annotation.
# The difference between this and 'add_gs' is that here, if this sequence and
# feature have an annotation already, it will be overwritten instead of appended to.
#
# Usage:
#        $stock->set_gs ($seq, $gsFeature, $line);
#        $stock->set_gs ($seq, $gsFeature, $line1, $line2, ...);
#        $stock->set_gs ($seq, $gsFeature, @lines);
#        $stock->set_gs ($seq, $gsFeature, "$line1\n$line2\n$...");
#
# Each list entry gets printed on its own line.
#
sub set_gs {
  my ($self, $seq, $feature, @text) = @_;
  unless (@text > 0) {
      carp "Missing parameters to 'set_gs' method; nothing will be done";
      return;
  }
  @text = map {$_ eq '' ? '' : split ("\n", $_, -1)} @text;
  $self->{gs}->{$feature}->{$seq} = [@text];
}

=head2 get_gs

    my $gs = $stock->get_gs (SEQNAME, FEATURE)

Concatenates and returns all lines beginning '#=GS SEQNAME FEATURE ...'

=cut

# Get a #=GS annotation.
#
sub get_gs {
  my ($self, $seq, $feature) = @_;
  unless (defined $feature) {
      carp "Missing parameters to 'get_gs' method; nothing will be done";
      return undef;
  }
  if (defined $self->{gs}->{$feature}->{$seq}) {
      return join ("\n", @{$self->{gs}->{$feature}->{$seq}});
  }
  else {
      #carp "Annotation for sequence $seq, feature $feature not found; returning empty string";
      return undef;
  }
}

=head2 drop_columns

    my $success = $stock->drop_columns (@columns)

Drops a set of columns from the alignment and the GC and GR annotations.  Note
that this may break annotations where the column annotations are not
independent (e.g. you might drop one base in an RNA base pair, but not the
other).  The caller is responsible for making sure the set of columns passed
in does not mess up the annoation.

Arguments: a list of scalar, zero-based column indices to drop.

Returns success of operation (1 for success, 0 for failure)

=cut

# Drops a set of columns from the alignment and the GC and GR annotations.  Note
# that this may break annotations where the column annotations are not
# independent (e.g. you might drop one base in an RNA base pair, but not the
# other).  The caller is responsible for making sure the set of columns passed
# in does not mess up the annoation.
#
# Arguments: a list of scalar, zero-based column indices to drop.
#
# Returns success of operation (1 for success, 0 for failure)
#
sub drop_columns {
  my ($self, @columns) = @_;

  unless ($self->is_flush()) {
    carp "File is not flush; cannot (should not) drop columns";
    return 0;
  }

  # drop from sequences and per-column annotations
  foreach my $seq_or_annot ($self->{seqdata}, $self->{gc})
    {
      # for each sequence or annotation string
      foreach my $key (keys %$seq_or_annot)
	{
	  # convert string to array form
	  my @data = split (//, $seq_or_annot->{$key});
	  foreach my $col (@columns)
	    {
	      if ( ($col < 0) or ($col >= @data) ) {
		carp "Column $col is out of range [1, ", scalar(@data),
		  "] in $key; can't drop, skipping";
	      }
	      else {
		# blank out the column that should get dropped
		$data[$col] = '';
	      }
	    }
	  # re-assemble back to string
	  $seq_or_annot->{$key} = join ('', @data);
	}
    }

  # drop from per-column/per-sequence annotations
  foreach my $feat_key (keys %{$self->{gr}})
    {
      foreach my $seq_key (keys %{$self->{gr}->{$feat_key}})
	{
	  my @data = split (//, $self->{gr}->{$feat_key}->{$seq_key});
	  foreach my $col (@columns)
	    {
	      if ( ($col < 0) or ($col >= @data) ) {
		carp "Column $col is out of range [1, ", scalar(@data),
		  "] in $feat_key, $seq_key; can't drop, skipping";
	      }
	      else {
		# blank out the column that should get dropped
		$data[$col] = '';
	      }
	    }
	  # re-assemble back to string
	  $self->{gr}->{$feat_key}->{$seq_key} = join ('', @data);
	}
    }

  return 1;
}

# parse line of Stockholm format file
# returns true if it finds a separator
sub parse_input_line {
    my ($self, $line) = @_;
    $line = $_ unless defined $line;

    # "#=GF [feature] [data]" line
    if ($line =~ /^\s*\#=GF\s+(\S+)\s+(\S.*)\s*$/) 
    { 
	my ($feature, $data) = ($1, $2);
	$self->add_gf ($feature, $data);  # preserve order of #=GF lines, for crazy context-sensitive Stockholm semantics
    }

    # "#=GC [feature] [data]" line
    elsif ($line =~ /^\s*\#=GC\s+(\S+)\s+(\S+)\s*$/) 
    { 
	my ($feature, $data) = ($1, $2, $3);
	$self->gc->{$feature} .= $data; 
    }

    # "#=GS [seqname] [feature] [data]" line
    elsif ($line =~ /^\s*\#=GS\s+(\S+)\s+(\S+)\s+(\S.*)\s*$/)
    { 
	my ($seqname, $feature, $data) = ($1, $2, $3);
	my $gs = $self->gs_($feature);
	$gs->{$seqname} = [] unless exists $gs->{$seqname};
	push @{$gs->{$seqname}}, $data;
    }

    # "#=GR [seqname] [feature] [data]" line
    elsif ($line =~ /^\s*\#=GR\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
    { 
	my ($seqname, $feature, $data) = ($1, $2, $3);
 	my $gr = $self->gr_($feature);
	$gr->{$seqname} = "" unless exists $gr->{$seqname};
	$gr->{$seqname} .= $data; 
    }

    # Unrecognised "#" line
    elsif ($line =~ /^\s*\#.*$/) {
#	warn "Comment line $line";
    }

    # Alignment separator: return true to indicate loop exit
    elsif ($line =~ /^\s*\/\/\s*$/) 
    { 
	return 1;
    }

    # Sequence line
    elsif ($line =~ /^\s*(\S+)\s*(\S*)\s*$/) {  # allows empty sequence lines
	my ($seqname, $data) = ($1, $2);
	unless (defined $self->seqdata->{$seqname}) {
	    $self->seqdata->{$seqname} = "";
	    push @{$self->seqname}, $seqname;  # preserve sequence order, for tidiness
	}
	$self->seqdata->{$seqname} .= $data;
    }

    elsif ($line =~ /\S/) 
    { 
	carp "Ignoring line: $_";
    }

    # This line wasn't a alignment separator: return false
    return 0;
}

=head2 empty

    my $isEmpty = $stock->empty

Returns true if the alignment is empty (i.e. there are no sequences).

=cut

# empty?
sub empty {
    my ($self) = @_;
    return @{$self->seqname} == 0;
}

=head2 subseq

    my $subseq = $stock->subseq (SEQNAME, $startPos)
    my $subseq = $stock->subseq (SEQNAME, $startPos, $length)

Returns a subsequence of the named alignment row.

Note that the co-ordinates are with respect to the alignment columns (starting at zero), NOT sequence positions.
That is, gaps are counted.

If $length is omitted, this call returns the subsequence from $startPos to the end of the row.

=cut

# Returns a subsequence.
#
# Arguments:
#   $name - name/identifier of sequence
#   $start - 0-based index of the starting nuc
#   $len - length of subsequence (if not specified, returns subseq from $start to the end)
#
sub subseq {
  my ($self, $name, $start, $len) = @_;

  my $seq = $self->{seqdata}->{$name};

  unless ($seq) {
    carp "No sequence, returning empty string for subsequence of '$name'\n";
    return '';
  }

  my $seqlength = length $seq;

  if ( ($start > $seqlength) or ($start < 0) ) {
    carp
      "Subsequence start nuc $start is outside of sequence bounds ",
	"[0, " . $seqlength-1 . "]; returning empty string for subsequence of $name\n";
    return '';
  }

  if (defined $len) {
    carp "Specified length $len is < 1... are you sure you want this?\n" if $len < 1;
    return substr ($seq, $start, $len);
  }
  else {
    return substr ($seq, $start);
  }
}

=head2 get_column

    my $column = $stock->get_column ($colNum)

Extracts the specified column and returns it as a string.

The column number $colNum uses zero-based indexing.

=cut

# Extracts the specified column and returns it as a string.
#
# Arguments:
#   $col - the column number
#
sub get_column {
  my ($self, $col_num) = @_;
  my $col_as_string;

  unless ($self->is_flush()) {
    carp "File is not flush; cannot (should not) get column";
    return 0;
  }

  my $align_length = $self->columns();

  if ( ($col_num >= $align_length) or ($col_num < 0) ) {
    carp
      "Column index $col_num is outside of alignment bounds ",
	"[1, $align_length]; returning empty string for column $col_num\n";
    return '';
  }
  else {
    foreach my $seq (@{$self->{seqname}}) {
      $col_as_string .= substr ($self->{seqdata}->{$seq}, $col_num, 1);
    }
    return $col_as_string;
  }
}

=head2 map_align_to_seq_coords

    my ($start, $end) = $stock->map_align_to_seq_coords ($colStart, $colEnd, $seqName);
    die "Failed" unless defined ($start);

    my ($start, $end) = $stock->map_align_to_seq_coords ($colStart, $colEnd, $seqName, $seqStart);
    die "Failed" unless defined ($start);

Given a range of column indices ($colStart, $colEnd), returns the coordinates
of sequence $seqName contained within that range, accounting for gaps.

Returns start and end coordinates as a 2-ple list, or returns the empty list
if there is no sequence in the given range (i.e. range contains all gaps in
$seqName).

You can provide $seqStart, which is the coordinate of the first non-gap
character in $seqName.  If you don't, the default is that the sequence starts
at 1.

=cut

# Given a range of column indices ($colStart, $colEnd), returns the coordinates
# of sequence $seqName contained within that range, accounting for gaps.
#
# Returns start and end coordinates as a 2-ple list, or returns the empty list
# if there is no sequence in the given range (i.e. range contains all gaps in
# $seqName).
#
# You can provide $seqStart, which is the coordinate of the first non-gap
# character in $seqName.  If you don't, the default is that the sequence starts
# at 1.
#
sub map_align_to_seq_coords {
  my ($self, $colStart, $colEnd, $seqName, $seqStart) = @_;
  $seqStart = 1 unless defined $seqStart;

  my $seq = $self->subseq ($seqName, $colStart, ($colEnd - $colStart + 1));

  if (my @seq = $seq =~ /[^$gapCharsRegexp]/g)
  {
    my $rangeStart;

    if ($colStart == 0) {
      # special case: range starts in the first column
      $rangeStart = $seqStart;
    }
    else {
      my $prefix = $self->subseq ($seqName, 0, $colStart);
      $rangeStart = $seqStart + length (join ('', ($prefix =~ /[^$gapCharsRegexp]/g)));
    }

    return ($rangeStart, ($rangeStart + @seq - 1));
  }
  else {
    # no sequence in desired range (i.e. all gaps)
    return ();
  }
}

=head2 is_gap

    my $isGap = $stock->is_gap (SEQNAME, $colNum)

Return true if a given (SEQNAME,column) co-ordinate is a gap.

=cut

# Return true if a given (rowname,column) co-ordinate is a gap
sub is_gap {
    my ($self, $row_name, $col) = @_;
    my $s = $self->seqdata->{$row_name};
    my $c = substr ($s, $col, 1);
    return $c =~ /[$gapCharsRegexp]/;
}

=head2 subalign

    my $substock = $stock->subalign ($start, $len)
    my $substock = $stock->subalign ($start, $len, $strand)

Returns a Stockholm object representing a sub-alignment of the current alignment,
starting at column $start (zero-based) and containing $len columns.

If $strand is supplied and is equal to '-', then the sub-alignment will be reverse-complemented.

=cut


# subalignment accessor
sub subalign {
    my ($self, $start, $len, $strand) = @_;
    $strand = '+' unless defined $strand;
    my $stock = $self->copy;

    # Sequence data
    foreach my $seqname (@{$stock->seqname}) {
	$stock->seqdata->{$seqname} = safe_substr ($stock->seqdata->{$seqname}, $start, $len);
	if ($strand eq '-') {
	    $stock->seqdata->{$seqname} = revcomp ($stock->seqdata->{$seqname});
	}
    }

    # GC lines
    foreach my $tag (keys %{$stock->gc}) {
	$stock->gc->{$tag} = safe_substr ($stock->gc->{$tag}, $start, $len);
	if ($strand eq '-') {
	    $stock->gc->{$tag} = reverse $stock->gc->{$tag};
	}
    }

    # GR lines
    foreach my $tag (keys %{$stock->gr}) {
	foreach my $seqname (keys %{$stock->gr->{$tag}}) {
	    $stock->gr->{$tag}->{$seqname} = safe_substr ($stock->gr->{$tag}->{$seqname}, $start, $len);
	    if ($strand eq '-') {
		$stock->gr->{$tag}->{$seqname} = reverse ($stock->gr->{$tag}->{$seqname});
	    }
	}
    }

    return $stock;
}

=head2 concatenate

    $stock->concatenate ($stock2)

Concatenates another alignment onto the end of this one.

=cut

# concatenate another alignment onto the end of this one
sub concatenate {
    my ($self, $stock) = @_;

    # get widths
    my $cols = $self->columns;
    my $catcols = $stock->columns;

    # Sequence data
    foreach my $seqname (@{$stock->seqname}) {
	if (exists $self->seqdata->{$seqname}) {
	    $self->seqdata->{$seqname} .= $stock->seqdata->{$seqname};
	} else {
	    push @{$self->seqname}, $seqname;
	    $self->seqdata->{$seqname} = "." x $cols . $stock->seqdata->{$seqname};
	}
    }
    foreach my $seqname (@{$self->seqname}) {
	unless (exists $stock->seqdata->{$seqname}) {
	    $self->seqdata->{$seqname} .= "." x $catcols;
	}
    }

    # GC lines
    foreach my $tag (keys %{$stock->gc}) {
	if (exists $self->gc->{$tag}) {
	    $self->gc->{$tag} .= $stock->gc->{$tag};
	} else {
	    $self->gc->{$tag} = "." x $cols . $stock->gc->{$tag};
	}
    }
    foreach my $tag (keys %{$self->gc}) {
	unless (exists $stock->gc->{$tag}) {
	    $self->gc->{$tag} .= "." x $catcols;
	}
    }

    # GR lines
    foreach my $tag (keys %{$stock->gr}) {
	if (exists $self->gr->{$tag}) {
	    foreach my $seqname (keys %{$stock->gr->{$tag}}) {
		if (exists $self->gr->{$tag}->{$seqname}) {
		    $self->gr->{$tag}->{$seqname} .= $stock->gr->{$tag}->{$seqname};
		} else {
		    $self->gr->{$tag}->{$seqname} = "." x $cols . $stock->gr->{$tag}->{$seqname};
		}
	    }
	} else {
	    $self->gr->{$tag} = {};
	    foreach my $seqname (keys %{$stock->gr->{$tag}}) {
		$self->gr->{$tag}->{$seqname} = "." x $cols . $stock->gr->{$tag}->{$seqname};
	    }
	}
    }
    foreach my $tag (keys %{$self->gr}) {
	foreach my $seqname (keys %{$self->gr->{$tag}}) {
	    unless (exists $stock->gr->{$tag} && exists $stock->gr->{$tag}->{$seqname}) {
		$self->gr->{$tag}->{$seqname} .= "." x $catcols;
	    }
	}
    }

    # GF and GS lines
    foreach my $tag (keys %{$stock->gf}) {
	$self->gf->{$tag} = [] unless exists $self->gf->{$tag};
	push @{$self->gf->{$tag}}, @{$stock->gf->{$tag}};
    }

    foreach my $tag (keys %{$stock->gs}) {
	$self->gs->{$tag} = {} unless exists $self->gs->{$tag};
	foreach my $seqname (keys %{$stock->gs->{$tag}}) {
	    $self->gs->{$tag}->{$seqname} = [] unless exists $self->gs->{$tag}->{$seqname};
	    push @{$self->gs->{$tag}->{$seqname}}, @{$stock->gs->{$tag}->{$seqname}};
	}
    }
}

# Reverse complement a sequence
sub revcomp {
  my ($arg0, $arg1) = @_;
  my $seq;

  if (defined ($arg1)) {  # must have been called externally (as method or package sub)
    $seq = $arg1;
  } else {
    $seq = $arg0;
  }

  $seq =~ tr/acgtuACGTU/tgcaaTGCAA/;
  $seq = reverse $seq;
  return $seq;
}

# Determine maximum sequence name length
sub maxNameLen {
    my ($self) = @_;
    return max(map(length($_), @{$self->seqname}));
}

=head2 is_flush

    my $flush = $stock->is_flush()

Returns true if the alignment is "flush", i.e. all sequence and annotation lines have the same length.

=cut

# Is the file flush? (i.e., are all sequences and annotations the same length?)
#
sub is_flush {
  my ($self) = @_;
  my $columns = $self->columns();

  # do stuff that has one key

  map {
    if (length($self->{seqdata}->{$_}) != $columns) {
      carp "Sequence $_ is not flush!";
      return 0;
    }
  } keys %{$self->{seqdata}};

  map {
    if (length($self->{gc}->{$_}) != $columns) {
      carp "#=GC annotation for feature $_ is not flush!";
      return 0;
    }
  } keys %{$self->{gc}};

  # do stuff that has two keys

  map {
    if (length($_) != $columns) {
      carp "#=GR annotation is not flush! ($_)";
      return 0;
    }
  } map { values %{$self->{gr}->{$_}} } keys %{$self->{gr}};

  return 1;  # if we got this far, everything must be flush
}

# Determine the list maximum value.
sub max  {
  my ($x, @y) = @_;
  foreach my $y (@y) 
  { 
    $x = $y if !defined($x) || (defined($y) && $y > $x)
  }
  return $x;
}

# Get/set the New Hampshire tree
sub NH {
    my ($self, $newval) = @_;
    if (defined $newval) {
	@{$self->gf_NH} = ($newval);
    }
    return join ("", @{$self->gf_NH});
}


=head2 add_row

    my $new_row_index = $stock->add_row ($seqname, $rowdata);

Adds a new row to a Stockholm alignment.

=cut

# add_row method
sub add_row {
    my ($self, $seqname, $rowdata) = @_;
    die "Attempt to add duplicate row" if exists $self->seqdata->{$seqname};
    push @{$self->seqname}, $seqname;
    $self->seqdata->{$seqname} = $rowdata;
    return @{$self->seqname} - 1;  # new row index
}

1;
