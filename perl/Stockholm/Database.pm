#!/usr/bin/perl -w

=head1 NAME

Stockholm::Database.pm

=head1 SYNOPSIS

Lightweight Perl module encapsulating a flatfile database of Stockholm multiple alignments.

For more detail on Stockholm format itself, see the following URL:

http://biowiki.org/StockholmFormat

=head1 GENERAL USAGE

A Stockholm::Database object is a blessed array reference.
(This is in contrast to most Perl objects, which are blessed hash references.)

The referenced array contains the constituent Stockholm objects.
So, if $db is a Stockholm::Database, then @$db is simply an array of Stockholm objects.

=head1 METHODS

=cut

package Stockholm::Database;

use strict;
use vars '@ISA';

use Stockholm;
use Carp;

=head2 new

    my $db = Stockholm::Database->new();

Creates an empty Stockholm::Database object.

=cut

# constructor
sub new {
    my ($class) = @_;
    my $self = [];
    bless ($self, ref($class) ? ref($class) : $class);
    return $self;
}

=head2 align_callback

    Stockholm::Database->align_callback ("RfamSeedAlignments.stock",
    				         sub
					 {
					     my $stock = shift;  # Stockholm object
					     # ...do something with the alignment...
					 })

Streams through alignments in a file, calling a callback function on each alignment as it is parsed.

This method can be used as an efficient alternative to reading the entire database into memory with the from_file method, for some applications.

=cut

# align_callback method
# streams through alignments in a file, calling a callback function as it parses each one
# can use in place of from_file method for some applications
sub align_callback {
    my ($class, $filename, $callback) = @_;

    my $self = $class->new;
    my $stock;

    local *FILE;
    local $_;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    while (<FILE>) {
	next unless /\S+/;  # skip blank lines
	$stock = Stockholm->new unless defined $stock;
	if ($stock->parse_input_line ($_)) {
	    &$callback ($stock);
	    $stock = undef;
	}
    }
    close FILE;
    &$callback ($stock) if defined $stock;

    return $self;
}

=head2 from_file

    my $db = Stockholm::Database->from_file ($filename)
    my $db = Stockholm::Database->from_file ($filename, 0)

Creates a Stockholm::Database object and populates it from a named Stockholm flatfile.

Reading a large database (e.g. PFAM) can take a while.
By default, a message is printed when each alignment is loaded.
The second form of the method suppresses these verbose messages.

=cut

# from_file method
sub from_file {
    my ($class, $filename, $verbose) = @_;
    $verbose = 1 unless defined $verbose;  # default

    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    my $self = $class->from_filehandle (\*FILE, $verbose);
    close FILE;

    return $self
}

=head2 from_filehandle

    my $db = Stockholm::Database->from_filehandle ($filehandle)
    my $db = Stockholm::Database->from_filehandle ($filehandle, 0)

Creates a Stockholm::Database object and populates it from a Stockholm flatfile, pointed to by a given filehandle.

Reading a large database (e.g. PFAM) can take a while.
By default, a message is printed when each alignment is loaded.
The second form of the method suppresses these verbose messages.

=cut

# from_filehandle method
sub from_filehandle {
    my ($class, $filehandle, $verbose) = @_;
    $verbose = 1 unless defined $verbose;  # default

    my $self = $class->new;
    my $stock;
    local $_;
    while (<$filehandle>) {
	next unless /\S+/;  # skip blank lines
	$stock = Stockholm->new unless defined $stock;
	if ($stock->parse_input_line ($_)) {
	    $self->add_alignment ($stock, $verbose);
	    $stock = undef;
	}
    }
    $self->add_alignment ($stock, $verbose) if defined $stock;

    return $self;
}


=head2 from_string

    $db->from_string ($string)
    $db->from_string ($string, 0)

Creates a Stockholm::Database object and populates it from a Stockholm-format string.

Reading a large database (e.g. PFAM) can take a while.
By default, a message is printed when each alignment is loaded.
The second form of the method suppresses these verbose messages.

=cut

# from_string
sub from_string {
    my ($class, $text, $verbose) = @_;
    my @text = map (split(/\n/), $text);

    my $self = $class->new;
    my $stock;

    foreach my $line (@text) {
	next unless $line =~ /\S+/;  # skip blank lines
	$stock = Stockholm->new unless defined $stock;
	if ($stock->parse_input_line ($line)) {
	    $self->add_alignment ($stock, $verbose);
	    $stock = undef;
	}
    }
    $self->add_alignment ($stock, $verbose) if defined $stock;

    return $self;
}

=head2 to_file

    $db->to_file ($filename)
    $db->to_file ($filename, $maxcols)

Saves a Stockholm::Database to a named file.

The second form of the method sets the maximum column width for the Stockholm::to_string method.

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

    $db->to_string()
    $db->to_string ($maxcols)

Renders a Stockholm::Database as a (Stockholm-format) string.

The second form of the method sets the maximum column width for the Stockholm::to_string method.

=cut

# to_string
sub to_string {
    my ($self, $maxcols) = @_;
    return join ("", map ($_->to_string ($maxcols), @$self));
}


=head2 add_alignment

    $db->add_alignment ($stock)
    $db->add_alignment ($stock, 0)

Adds a Stockholm alignment to a Stockholm::Database.

By default, a message is printed when the alignment is added.
The second form of the method suppresses these verbose messages.

=cut

# add_alignment method
sub add_alignment {
    my ($self, $stock, $verbose) = @_;
    push @$self, $stock;
    my $id = $stock->get_gf ('ID');
    warn "...loaded Stockholm alignment", defined $id ? " '$id'" : '', " (", $stock->sequences, " sequences)\n" if $verbose;
    return $self;
}

# end of package

1;
