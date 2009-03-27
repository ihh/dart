
=head1 NAME

DartSexpr.pm

=head1 SYNOPSIS

Perl module for working with simple, DART-style S-expressions.

=head1 EXAMPLES

=head2 Example #1

  use DartSexpr;
  my $e = DartSexpr->new(["hello","there"],["world",["beauty",[qw(trees roses skies clouds)]]]);
  print
    $e->tag_value ("hello"),
    "\n",
    $e->hello->value("beautiful"),
    "\n",
    $e->to_string,
    "\n";

=head2 Example #2

  use DartSexpr;
  print
    DartSexpr->from_string("(man big) (woman kind) (baby cute)")->woman->value,
    "\n";

=head2 Example #3

  use DartSexpr;
  my $filename = "$ENV{DARTDIR}/data/codon.eg";
  print DartSexpr->from_file($filename)->alphabet->to_string,"\n";

=head1 GENERAL USAGE

A DartSexpr object is a blessed array reference.
(This is in contrast to most Perl objects, which are blessed hash references.)

The structure is recursive:
any of the elements of the array can either be DartSexpr objects themselves,
or they can be strings (i.e. "atoms" in Lisp parlance).

=head1 METHODS

=cut

package DartSexpr;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw (new clone add add_child set
	      grep grep_child map_child
	      find_all find find_or_add
	      has_tag tag value values tag_value tag_values
	      from_string to_string from_file to_file
	      AUTOLOAD recursive_bless_array_refs);
@EXPORT_OK = @EXPORT;

use strict;
use vars '@ISA';

use Carp;

=head2 new

    my $sexpr = DartSexpr->new();
    my $sexpr = DartSexpr->new(['a','b',['c','d']]);

Creates a new S-expression.
Array references supplied as args are auto-blessed into S-expressions too.

=cut

# constructor: creates a new S-expression.
# Array references supplied as args are auto-blessed into S-expressions too.
sub new {
    my ($class, @data) = @_;
    my $self = [@data];
    $class = ref($class) if ref($class);
    $class->recursive_bless_array_refs ($self);
    return $self;
}

=head2 clone

    my $sexpr2 = $sexpr1->clone();

Deep-copy constructor.

=cut

# deep copy constructor
sub clone {
    my ($self) = @_;
    my $clone = ref($self)->new;
    $clone->add (map (ref($_) ? $_->clone : $_, @$self));
    return $clone;
}

# recursive_bless_array_refs method
sub recursive_bless_array_refs {
    my ($class, @array) = @_;

    # bless array refs
    @array = grep (ref($_) eq 'ARRAY', @array);
    foreach my $array (@array) { bless $array, $class }

    # call this subroutine on child arrays
    @array = map (@$_, @array);
    $class->recursive_bless_array_refs (@array) if @array;
}

=head1 ACCESSOR METHODS

=head2 Autoloaded accessor

    $sexpr->TAG
    $sexpr->TAG($newValue)

This is a shorthand that is [almost] equivalent to the find_or_add method, i.e.

    $sexpr->find_or_add(TAG)
    $sexpr->find_or_add(TAG,$newValue)

The "[almost]" qualifier in the above refers to the following awkward hack:
if TAG contains any underscores, they are transliterated to hyphens before the call to find_or_add.

For example, this...

    $sexpr->update_policy ('irrev')

...is exactly equivalent to this...

    $sexpr->find_or_add ('update-policy', 'irrev')

=head2 add

    $sexpr->add(@data);

Adds an element (or elements) to the top-level list in the S-expression.
The added elements can include references to arrays, which are then recursively blessed into DartSexpr objects.

=cut

# add method: adds an element (or elements) to the list
sub add {
    my ($self, @data) = @_;
    @data = map (ref($_) ? $_ : split(/\s+/,$_), @data);
    ref($self)->recursive_bless_array_refs (@data);
    push @$self, @data;
    return $self;
}

=head2 add_child

    my $child = $sexpr->add_child()
    my $child = $sexpr->add_child(@data)

Creates a new child element and adds it to the list.
The array @data is added to the child.

=cut

# add_child method: creates a new child element and adds it to the list
sub add_child {
    my ($self, @data) = @_;
    my $child = DartSexpr->new;
    push @$self, $child;
    $child->add (@data);
    return $child;
}

=head2 set

    my $child_sexpr = $sexpr->set (NAME, @data)

Finds a named child node (using the find_or_add method),
then changes that node\'s value(s) to @data.

Example:

    use DartSexpr
    $s = DartSexpr->new ([a,b], [c,d], [e,f]);
    print $s->to_string(),"\n";

(a b) (c d) (e f)

    $s->set ("c", qw(x y));
    print $s->to_string(),"\n";

(a b) (x y) (e f)

=cut

# set method: changes a named node
sub set {
    my ($self, $name, @data) = @_;
    croak "name undefined" unless defined $name;
    @data = map (ref($_) ? $_ : split(/\s+/,$_), @data);
    ref($self)->recursive_bless_array_refs (@data);
    my $sexpr = $self->find_or_add ($name);
    @$sexpr = @data;
    return $sexpr;
}

# By-name accessors.

=head2 grep

    my @match_sexpr = $sexpr->grep (sub {
	my $descendant_sexpr = shift;
	...
	# return true or false
    });

Returns all elements matching a predicate, recursively traversing the entire S-expression tree rooted at this node.

=cut

# grep method: returns all elements matching a predicate
sub grep {
    my ($self, $predicate) = @_;
    croak "bad predicate" unless ref($predicate) eq 'CODE';
    return grep (&$predicate($_), @$self);
}

=head2 grep_child

    my @match_sexpr = $sexpr->grep_child (sub {
	my $child_sexpr = shift;
	...
	# return true or false
    });

Returns all elements matching a predicate, visiting only the immediate children of this node.

=cut

# grep_child method: returns all children matching a predicate
sub grep_child {
    my ($self, $predicate) = @_;
    croak "bad predicate" unless ref($predicate) eq 'CODE';
    return grep (ref($_) && &$predicate($_), @$self);
}

=head2 map_child

    my @mapped_sexpr = $sexpr->map_child (sub {
	my $child_sexpr = shift;
	...
	# return some function of $child_sexpr
    });

Returns the application of a given function to the list of immediate children of this node.

=cut

# map_child method: returns application of given function to all children
sub map_child {
    my ($self, $func) = @_;
    croak "bad map function" unless ref($func) eq 'CODE';
    return map (ref($_) ? &$func($_) : (), @$self);
}

=head2 find_all

    my @found_sexpr = $sexpr->find_all ($tag);

Searches the immediate children of this node, looking for nodes (a) that are lists and (b) whose first element is the atom $tag.

=cut

# find_all method: returns all children whose first element matches tag
sub find_all {
    my ($self, $tag) = @_;
    croak "tag undefined" unless defined $tag;
    return $self->grep_child (sub {
	my ($child) = @_;
	return $child->has_tag && $child->tag eq $tag;
    });
}

=head2 find

    my $found_sexpr = $sexpr->find ($tag);

Searches the immediate children of this node, looking for nodes (a) that are lists and (b) whose first element is the atom $tag.

Returns exactly one such node, or generates an error.

=cut

# find method: returns unique child whose first element matches name, or dies
sub find {
    my ($self, $name) = @_;
    croak "name undefined" unless defined $name;
    my @found = $self->find_all ($name);
    croak "More than one child named $name" if @found > 1;
    croak "In S-expression '", $self->to_string, "':\nNo child named $name" if @found < 1;
    return $found[0];
}

=head2 find_or_add

    my $found_sexpr = $sexpr->find_or_add ($tag)
    my $found_sexpr = $sexpr->find_or_add ($tag, $newValue)
    my $found_sexpr = $sexpr->find_or_add ($tag, @newValues)

Searches the immediate children of this node, looking for nodes (a) that are lists and (b) whose first element is the atom $tag.

If no such node exists, one is added and returned.

If more than one such node exists, the method generates an error.

The latter forms set the new "value" or "values" of the node,
i.e. the elements of the list that follow the "tag".

=cut

# find_or_add method: finds (or creates) a child whose first element matches name.
# Dies if more than one such child exists.
# Returns the child element.
sub find_or_add {
    my ($self, $name, @data) = @_;
    my @found = $self->find_all ($name);
    croak "More than one child named $name" if @found > 1;
    if (@found < 1) {
	carp "...auto-adding '$name' to S-expression";
	return $self->add_child ($name, @data);
    }
    return $found[0];
}

# AUTOLOAD delegates to find_or_add or set, depending on whether arguments are supplied
sub AUTOLOAD {
    my ($self, @args) = @_;
    my $sub = our $AUTOLOAD;
    $sub =~ s/.*:://;  # strip off module path

    # check for DESTROY
    return if $sub eq "DESTROY";

    # pass to find or set
    $sub =~ tr/_/-/;  # replace underscores with hyphens for cute S-expressions
    return @args
	? $self->set ($sub, $sub, @args)  # prepend the tagname (can't use AUTOLOAD to change this)
	: $self->find_or_add ($sub);
}

# Accessors for children of the form (tag value) or (tag value1 value2...)

=head2 has_tag

    my $has_tag = $sexpr->has_tag()

Returns true if the node is a list and its first element is an atom.

=cut

# has_tag
sub has_tag {
    my ($self) = @_;
    return @$self && !ref($self->[0]);
}

=head2 tag

    my $tag = $sexpr->tag()

Returns the first element of the list, which must be an atom.

If the list is empty, or the first element is not an atom, it generates an error.

=cut

# tag
sub tag {
    my ($self) = @_;
    croak "S-expression is empty" unless @$self;
    croak "Tag undefined: first child is not a string" if ref($self->[0]);
    return $self->[0];
}

=head2 values

    my @values = $sexpr->values()
    $sexpr->values (@new_values)

Get/set the "values" of this node, i.e. all elements except the first.

If the node does not have a "tag"
(i.e. the list is empty, or the first element is not an atom), it generates an error.

=cut

# values accessor, for tags with (potentially) multiple values
sub values {
    my ($self, @values) = @_;
    croak "Node is not a (tag value[s]) pair: there's no tag" if @$self < 1;
    croak "Node (", $self->tag, " ...) is not a (tag value[s]) pair: there's no value" if @$self < 2 && !@values;
    @$self = ($self->tag, @values) if @values;
    return @$self[1..@$self-1];
}

=head2 value

    my $value = $sexpr->value()
    $sexpr->value ($new_value)

Get/set the "value" of this node, i.e. the second element of the list.

If the node does not have a "tag"
(i.e. the list is empty, or the first element is not an atom), it generates an error.

=cut

# value accessor, for tags with single values
sub value {
    my ($self, $value) = @_;
    my @values = $self->values (defined($value) ? ($value) : ());
    croak "Node is not a (tag value) pair: there's more than one value" if @values > 1;
    return $values[0];
}

=head2 tag_values

    my @tag_values = $sexpr->tag_values()
    $sexpr->tag_values ($new_tag, @new_values)

Get/set both the "tag" and the "values" of this node.

=cut

# tag_values accessor, for tags with (potentially) multiple values
sub tag_values {
    my ($self, $tag, @values) = @_;
    croak "tag undefined" unless defined $tag;
    my $sexpr = @values
	? $self->find_or_add ($tag, @values)
	: $self->find($tag);
    return $sexpr->values (@values);
}

=head2 tag_value

    my $tag_value = $sexpr->tag_value()
    $sexpr->tag_value ($new_tag, $new_value)

Get/set both the "tag" and the "value" of this node.

=cut

# tag_value accessor, for tags with single values
sub tag_value {
    my ($self, $tag, $value) = @_;
    croak "tag undefined" unless defined $tag;
    my $sexpr = defined($value)
	? $self->find_or_add ($tag, $value)
	: $self->find($tag);
    return $sexpr->value ($value);
}

# Input/output methods.

=head1 INPUT/OUTPUT METHODS

=head2 to_string

    print $sexpr->to_string()

Prints the S-expression as a string.

=cut

# to_string method
sub to_string {
    my ($self) = @_;
    return join (" ", map (ref($_) ? "(" . $_->to_string . ")" : $_, @$self));
}

=head2 from_string

    my $sexpr = DartSexpr->from_string($string)
    my $sexpr = DartSexpr->from_string(@strings)

Creates a new S-expression and parses its contents from a string.

=cut

# from_string method
sub from_string {
    my ($class, @strings) = @_;
    my $self = $class->new;

    grep (s/;;.*$//, @strings);
    my $string = join ("", @strings);
    $string =~ s/[\s\r]+/ /g;

    my @split;
    my $l = undef;
    my $nest = 0;

    for (my $i = 0; $i < length($string); ++$i) {
	my $c = substr ($string, $i, 1);

	if ($c eq "(") {  # opening bracket
	    if ($nest++ == 0) {
		$self->add (substr ($string, $l, $i - $l)) if defined $l;
		$l = $i;
	    }

	} elsif ($c eq ")") {  # closing bracket
	    croak "$string\nToo many )'s" if $nest <= 0;
	    if (--$nest == 0) {
		$self->add (DartSexpr->from_string (substr ($string, $l + 1, $i - 1 - $l)));
		$l = undef;
	    }

	} elsif ($nest == 0) {
	    if ($c =~ /\S/) {  # non-whitespace
		$l = $i unless defined $l;

	    } else {  # whitespace
		$self->add (substr ($string, $l, $i - $l)) if defined $l;
		$l = undef;
	    }
	}
    }
    croak "$string\nToo many ('s" if $nest;
    $self->add (substr ($string, $l, length($string) - $l)) if defined $l;

    return $self;
}

=head2 from_file

    my $sexpr = DartSexpr->from_file($filename)

Creates a new S-expression and parses its contents from a file.

=cut

# from_file method
sub from_file {
    my ($class, $filename) = @_;

    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open '$filename': $!";
    my @file = <FILE>;
    close FILE;

    print STDERR "...parsing S-expression file '$filename'... ";
    my $self = $class->from_string (@file);
    print STDERR "done\n";

    return $self;
}

=head2 from_filehandle

    my $sexpr = DartSexpr->from_filehandle($filehandle)

Creates a new S-expression and parses its contents from a filehandle.

=cut

sub from_filehandle {
    my ($class, $filehandle) = @_;

    my @file = <$filehandle>;
    close $filehandle;

    print STDERR "...parsing S-expression file... ";
    my $self = $class->from_string (@file);
    print STDERR "done\n";

    return $self;
}

=head2 to_file

    $sexpr->to_file ($filename)

Prints the S-expression to a named file.

=cut


# to_file method
sub to_file {
    my ($self, $filename) = @_;

    local *FILE;
    open FILE, ">$filename" or croak "Couldn't open '$filename' for writing: $!";
    print FILE $self->to_string;
    close FILE or croak "Couldn't close '$filename': $!";;

    return $filename;
}

# End of package
1;
