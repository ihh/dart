package Microarray::Data;

use Carp;
use strict qw(vars subs refs);

1;

# constructor

sub new {
    my ($class, $filename) = @_;

    my $self = { 'row_data' => [],
		 'column_name' => undef,
		 'row_name' => [],
		 'column_index' => {},
		 'row_index' => {},
		 '_got_bounds' => 0,
		 'min' => undef,
		 'max' => undef };
    bless $self, $class;

    local *FILE;
    open FILE, "<$filename" or confess $!;

    my $col_name = <FILE>;
    chomp $col_name;
    my ($dummy, @col_name) = split /\t/, $col_name;
    $self->{'column_name'} = \@col_name;
    for (my $c = 0; $c < @col_name; ++$c) {
	my $cl = $col_name[$c];
	confess "Duplicate column '$cl'" if exists $self->column_index->{$cl};
	$self->column_index->{$cl} = $c;
    }

    while (<FILE>) {
	chomp;
	my ($row_name, @row_data) = split /\t/;
	$self->add_row ($row_name, @row_data);
    }

    close FILE;

    $self;
}

# deep copy constructor

sub copy {
    my ($orig) = @_;
    my $self = { 'row_data' => [map ([@$_], @{$orig->row_data})],
		 'column_name' => [@{$orig->column_name}],
		 'row_name' => [@{$orig->row_name}],
		 'column_index' => {%{$orig->column_index}},
		 'row_index' => {%{$orig->row_index}},
		 '_got_bounds' => $orig->_got_bounds,
		 'min' => $orig->{'min'},
		 'max' => $orig->{'max'} };
    bless $self, ref $orig;
    $self;
}

# accessors

sub row_data {
    my ($self) = @_;
    $self->{'row_data'};
}

sub column_name {
    my ($self) = @_;
    $self->{'column_name'};
}

sub row_name {
    my ($self) = @_;
    $self->{'row_name'};
}

sub row_index {
    my ($self) = @_;
    $self->{'row_index'};
}

sub column_index {
    my ($self) = @_;
    $self->{'column_index'};
}

sub _got_bounds {
    my ($self, $_got_bounds) = @_;
    if (defined $_got_bounds) { return $self->{'_got_bounds'} = $_got_bounds }
    $self->{'_got_bounds'};
}

sub min {
    my ($self) = @_;
    $self->_get_bounds;
    $self->{'min'};
}

sub max {
    my ($self) = @_;
    $self->_get_bounds;
    $self->{'max'};
}

# output method

sub tab_delimit {
    my ($self) = @_;
    my $text = "\t" . join ("\t", @{$self->column_name}) . "\n";
    for (my $row = 0; $row < $self->rows; ++$row) {
	$text .= $self->row_name->[$row] . "\t" . join ("\t", map (defined() ? $_ : "", @{$self->row_data->[$row]})) . "\n";
    }
    $text;
}

# helpers

sub columns {
    my ($self) = @_;
    @{$self->column_name} + 0;
}

sub rows {
    my ($self) = @_;
    @{$self->row_name} + 0;
}

sub clear_rows {
    my ($self) = @_;
    $self->{'row_name'} = [];
    $self->{'row_index'} = {};
    $self->{'row_data'} = [];
    $self->_got_bounds (0);
}

sub add_row {
    my ($self, $row_name, @row_data) = @_;

    warn "Row '$row_name' has ", @row_data+0, " columns (should have ", $self->columns, ")" unless @row_data == $self->columns;
    @row_data = @row_data[0..$self->columns-1];
    foreach (@row_data) { $_ = undef if defined && !/\d/ }

    push @{$self->row_name}, $row_name;
    confess "Duplicate row '$row_name'" if exists $self->row_index->{$row_name};
    push @{$self->row_data}, \@row_data;
    $self->row_index->{$row_name} = $self->rows - 1;
    $self->_got_bounds (0);
}

sub get_row_by_name {
    my ($self, $row_name) = @_;
    my $row_index = $self->row_index->{$row_name};
    confess "Row '$row_name' not found" unless defined $row_index;
    $self->row_data->[$row_index];
}

# get bounds

sub _get_bounds {
    my ($self) = @_;
    if (!$self->_got_bounds) {

	my ($min, $max);

	foreach my $row_data (@{$self->row_data}) {
	    foreach my $x (@$row_data) {
		$min = $x if !defined ($min) || (defined ($x) && $x < $min);
		$max = $x if !defined ($max) || (defined ($x) && $x > $max);
	    }
	}

	$self->{'min'} = $min;
	$self->{'max'} = $max;

	$self->_got_bounds (1);
    }
}

# get cumulative bin positions
# on entry, $n_bins = number of bins
# returns an array @x such that a proportion ($i/$n_bins) of entries are <= $x[$i]

sub get_bins {
    my ($self, $n_bins, $keep_symmetric) = @_;

    my @sorted_neg = sort { $a <=> $b } map (grep (defined() && $_ < 0, @{$self->row_data->[$_]}), 0..$self->rows-1);
    my @sorted_pos = sort { $a <=> $b } map (grep (defined() && $_ >= 0, @{$self->row_data->[$_]}), 0..$self->rows-1);

    if (defined ($keep_symmetric) && $keep_symmetric) {
	if (-$sorted_neg[@sorted_neg-1] > $sorted_pos[@sorted_pos-1]) {
	    @sorted_pos = map (-$_, @sorted_neg);
	} else {
	    @sorted_neg = map (-$_, @sorted_pos);
	}
    }

    my @bin;
    for (my $bin = 0; $bin < ($n_bins-1)/2; ++$bin) {
	push @bin, $sorted_neg[@sorted_neg * $bin / ($n_bins/2)];
    }
    for (my $bin = 0; $bin < $n_bins/2; ++$bin) {
	push @bin, $sorted_pos[(@sorted_pos - 1) * $bin / ($n_bins/2)];
    }
    return @bin;
}

sub get_bins_ignore_sign {
    my ($self, $n_bins) = @_;

    my @sorted = sort { $a <=> $b } map (grep (defined(), @{$self->row_data->[$_]}), 0..$self->rows-1);
    my @bin;
    for (my $bin = 0; $bin < $n_bins; ++$bin) {
	push @bin, $sorted[@sorted * $bin / $n_bins];
    }
    return @bin;
}

# make a kimono covariance matrix from a supplied covariance function

sub covariance_matrix {
    my ($self, $cov_fn) = @_;
    if (!defined $cov_fn) { $cov_fn = sub { shift() eq shift() ? 1 : 0 } }  # default covariance matrix is identity
    my $cn = $self->column_name;
    my $text = "\t" . join ("\t", @$cn) . "\n";
    for (my $row = 0; $row < @$cn; ++$row) {
	$text .= $$cn[$row];
	for (my $col = 0; $col <= $row; ++$col) {
	    $text .= "\t" . &$cov_fn ($$cn[$row], $$cn[$col]);
	}
	$text .= "\n";
    }
    $text;
}

# offset
sub offset {
    my ($self) = @_;
    my $min = $self->min;
    my $max = $self->max;
    foreach my $row_data (@{$self->row_data}) {
	foreach my $entry (@$row_data) {
	    if (defined $entry) {
		$entry -= $min + ($max-$min)/2;
	    }
	}
    }
    $self->{'_got_bounds'} = 0;
}
