
=head1 NAME

PhyloGram::Chain.pm

=head1 SYNOPSIS

Perl module for working with Markov chain substitution models in xrate-format phylo-grammar files.

The module extends the DartSexpr.pm module, and so inherits all of the methods from that class.

For a description of the format itself, see the following URL:

http://biowiki.org/XrateFormat

=head1 METHODS

=cut

# Perl module for working with xgram format Markov chains

package PhyloGram::Chain;

use DartSexpr;

@ISA = qw(DartSexpr);

use strict;
use vars '@ISA';

use Carp;

=head2 new_chain

    my $chain = PhyloGram::Chain->new_chain (@pseudoterminals)

Creates a new PhyloGram::Chain object with the specified pseudoterminals.

The 'update-policy' is initialized to 'irrev'.

=cut

# new_chain method
# all rates & probabilities zero
sub new_chain {
    my ($class, @pseudoterm) = @_;
    my $self = DartSexpr->new ('chain',
			       ['update-policy' => 'irrev'],
			       ['terminal' => \@pseudoterm],
			       );
    bless $self, $class;
    return $self;
}

=head2 from_xrate

    my $chain = PhyloGram::Chain->from_xrate ($filename, $pseudoterminal)

Creates a new PhyloGram::Chain object and populates it from an xrate-format file,
 by scanning the file for a chain declaration containing the specified pseudoterminal.

=cut

# from_xrate method: creates a new, single-pseudoterminal chain from an xrate-format file
sub from_xrate {
    my ($class, $file, $pseudoterm) = @_;
    $pseudoterm = "default_xrate_pseudoterm" unless defined $pseudoterm;
    my $self = $class->new_chain ($pseudoterm);

    local *XRATE;
    local $_;
    open XRATE, "<$file" or croak "Couldn't open $file: $!";

    $_ = <XRATE>;  # CA line
    my ($c, $a) = split;
    croak "xrate file contains hidden classes, not allowed in chain" if $c != 1;

    $_ = <XRATE>;  # alphabet
    my @tok = split;

    $_ = <XRATE>;  # eqm freqs
    my @eqm = split;
    for (my $i = 0; $i < @tok; ++$i) {
	$self->initial ($tok[$i], $eqm[$i]);
    }

    # rate matrix
    my $mutate_hash = $self->mutate_hash;
    for (my $i = 0; $i < @tok; ++$i) {
	$_ = <XRATE>;
	my @rate = split;
	for (my $j = 0; $j < @tok; ++$j) {
	    if ($j != $i) {
		$self->mutate (@tok[$i,$j], $rate[$j], $mutate_hash);
	    }
	}
    }

    close XRATE;
    return $self;
}

=head2 initial

    my $initSexpr = $chain->initial ($state)
    my $initSexpr $chain->initial ($state, $newProb)

The get/set accessor for the initial probability for a state.
The latter form sets the new probability, creating the S-expression if it doesn\'t already exist.

For chains with more than one pseudoterminal, the $state can be a whitespace-separated string of alphabet tokens, or it can be a reference to an array of alphabet tokens.

In both forms of the method, the return value is the S-expression containing the initial probability definition, i.e.

    (initial (state ...) (prob ...))

If (as is often the case) what you want is the actual probability, use something like the following (see DartSexpr.pm documentation for explanations of how this works):

    $initSexpr->prob->value

=cut

# accessor for the initial probability for a state
# returns the S-expression (initial (state ...) (prob ...))
#  -- to get actual probability, use $e->prob->value;
# if $newProb is defined, then the S-expression will be created (if it doesn't already exist)
sub initial {
    my ($self, $state, $newProb) = @_;
    $state = [split(//,$state)] unless ref($state);

    my @sexpr = $self->grep_child (sub {
	my ($child) = @_;
	return 0 unless $child->has_tag && $child->tag eq 'initial';
	my $child_state = $child->state->value;
	return @$state == @$child_state && !grep ($$state[$_] ne $$child_state[$_], 0..@$state-1);
    });

    croak "More than one matching S-expression (initial (state (@$state)))" if @sexpr > 1;
    if (defined $newProb) {
	if (@sexpr) {
	    $sexpr[0]->prob($newProb);
	} else {
	    my $new_sexpr = DartSexpr->new ('initial',
					    ['state' => $state],
					    ['prob' => $newProb]);
	    $self->add ($new_sexpr);
	    push @sexpr, $new_sexpr;
	}
    }

    return $sexpr[0];
}

=head2 mutate

    my $mutSexpr = $chain->mutate ($srcState, $destState)
    my $mutSexpr = $chain->mutate ($srcState, $destState, $mutateHashRef)

    my $mutSexpr = $chain->mutate ($srcState, $destState, $newRate)
    my $mutSexpr = $chain->mutate ($srcState, $destState, $newRate, $mutateHashRef)

The get/set accessor for the mutation rate between two states.
The latter forms set the new rate, creating the S-expression if it doesn\'t already exist.

For chains with more than one pseudoterminal, the $srcState and $destState can be whitespace-separated strings of alphabet tokens, or they can be references to arrays of alphabet tokens.

If one of the forms with $mutateHashRef is used, where $mutateHashRef is a hash-reference returned by a call to the C<mutate_hash> method,
then the hash will be used to speed up element lookup.
This is HIGHLY RECOMMENDED especially for large chains, which can otherwise be VERY SLOW.

In all forms of the method, the return value is the S-expression containing the mutation rate definition, i.e.

    (mutate (from ...) (to ...) (rate ...))

If (as is often the case) what you want is the actual probability, use something like the following (see DartSexpr.pm documentation for explanations of how this works):

    $mutSexpr->rate->value

=cut

# accessor for the mutation rate between two states
# returns the S-expression (mutate (from ...) (to ...) (rate ...))
#  -- to get actual rate, use $e->rate->value;
# if $newRate is defined, then the S-expression will be created (if it doesn't already exist)
# if $mutate_hash is correctly defined (from a call to mutate_hash, below) then a hash will be used to speed up element access
# NB newRate can be omitted.
sub mutate {
    my ($self, $from, $to, $newRate, $mutate_hash) = @_;
    ($newRate, $mutate_hash) = (undef, $newRate) if ref($newRate);   # perl hackiness so caller can omit newRate
    $from = [split(//,$from)] unless ref($from);
    $to = [split(//,$to)] unless ref($to);
    croak "from & to have different lengths" unless @$from eq @$to;
# Commented out the following test because xrate chains now allow mutations-to-self... IH, 12/8/2008
#    croak "from==to" unless grep ($$from[$_] ne $$to[$_], 0..@$from-1);
    my $mutate_hash_key = join ("", @$from, ' ', @$to);

    my @sexpr;
    if (defined $mutate_hash) {
	@sexpr = exists ($$mutate_hash{$mutate_hash_key}) ? ($$mutate_hash{$mutate_hash_key}) : ();
    } else {
	@sexpr = $self->grep_child (sub {
	    my ($child) = @_;
	    return 0 unless $child->has_tag && $child->tag eq 'mutate';
	    my $child_from = $child->from->value;
	    my $child_to = $child->to->value;
	    return 0 unless @$from == @$child_from && !grep ($$from[$_] ne $$child_from[$_], 0..@$from-1);
	    return 0 unless @$to == @$child_to && !grep ($$to[$_] ne $$child_to[$_], 0..@$to-1);
	    return 1;
	});
    }

    croak "More than one matching S-expression (mutate (from (@$from)) (to (@$to)))" if @sexpr > 1;
    if (defined $newRate) {
	 if (@sexpr) {
	     $sexpr[0]->rate($newRate);
	 } else {
	     my $new_sexpr = DartSexpr->new ('mutate',
					     ['from' => $from],
					     ['to' => $to],
					     ['rate' => $newRate]);
	     $self->add ($new_sexpr);
	     push @sexpr, $new_sexpr;
	     $$mutate_hash{$mutate_hash_key} = $new_sexpr if defined $mutate_hash;
	 }
    }
    return $sexpr[0];
}

=head2 mutate_hash

    my $mutateHashRef = $chain->mutate_hash()

Returns a reference to a hashtable of all mutations in the chain.
This can be used to greatly speed up subsequent accesses of individual mutation rates via the C<mutate> method.

The hash is of little direct use to the caller, and is typically only passed to the C<mutate> method.

=cut

# method to return reference to a hashtable of all mutations in the chain
sub mutate_hash {
    my ($self) = @_;
    my %mutate;
    foreach my $mutate ($self->find_all ('mutate')) {
	$mutate{join("",@{$mutate->from->value},' ',@{$mutate->to->value})} = $mutate;
    }
    return \%mutate;
}

=head2 inc_initial

    $chain->inc_initial ($otherChain)
    $chain->inc_initial ($otherChain, $weight)

Adds (weighted) initial probabilities from another chain to this chain.

=cut

# method to add weighted initial probabilities from another chain
sub inc_initial {
    my ($self, $chain, $weight) = @_;
    $weight = 1 unless defined $weight;
    foreach my $chain_init ($chain->find_all ('initial')) {
	my $self_init = $self->initial ($chain_init->state->value);
	unless (defined $self_init) {
	    $self_init = $self->initial ($chain_init->state->value, 0);
	}
	$self_init->prob ($self_init->prob->value + $weight * $chain_init->prob->value);
    }
}

=head2 inc_rates

    $chain->inc_rates ($otherChain)
    $chain->inc_rates ($otherChain, $weight)

Adds (weighted) mutation rates from another chain to this chain.

=cut

# method to add weighted rates from another chain
sub inc_rates {
    my ($self, $chain, $weight) = @_;
    $weight = 1 unless defined $weight;
    my $mutate_hash = $self->mutate_hash;
    foreach my $chain_mutate ($chain->find_all ('mutate')) {
	my $self_mutate = $self->mutate ($chain_mutate->from->value, $chain_mutate->to->value, $mutate_hash);
	unless (defined $self_mutate) {
	    $self_mutate = $self->mutate ($chain_mutate->from->value, $chain_mutate->to->value, 0, $mutate_hash);
	}
	$self_mutate->rate ($self_mutate->rate->value + $weight * $chain_mutate->rate->value);
    }
}

=head2 scale_rates

    $chain->scale_rates ($scale_factor)

Scales all the rates in the chain by a constant factor.

=cut

# method to scale rates
sub scale_rates {
    my ($self, $scale_factor) = @_;
    my $mutate_hash = $self->mutate_hash;
    foreach my $self_mutate ($self->find_all ('mutate')) {
	$self_mutate->rate ($self_mutate->rate->value * $scale_factor);
    }
}

# end
1;
