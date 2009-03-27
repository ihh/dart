# Perl module for working with xgram grammar files in DNA

package PhyloGram::Dna;

use DartSexpr;
use PhyloGram;
use PhyloGram::Chain;

@ISA = qw(PhyloGram);

use strict;
use vars '@ISA';

use Carp;

# constructor
sub new {
    my ($class) = @_;

    # create S-expression
    my $self = PhyloGram->new_default_grammar;
    bless $self, $class;

    # add DNA alphabet
    $self->add (DartSexpr->from_string (<<END));
alphabet
(name DNA)
(token (a c g t))
(complement (t g c a))
(extend (to n) (from a) (from c) (from g) (from t))
(extend (to x) (from a) (from c) (from g) (from t))
(extend (to u) (from t))
(extend (to r) (from a) (from g))
(extend (to y) (from c) (from t))
(extend (to m) (from a) (from c))
(extend (to k) (from g) (from t))
(extend (to s) (from c) (from g))
(extend (to w) (from a) (from t))
(extend (to h) (from a) (from c) (from t))
(extend (to b) (from c) (from g) (from t))
(extend (to v) (from a) (from c) (from g))
(extend (to d) (from a) (from g) (from t))
(wildcard *)
END

    # return
    return $self;
}

# end of package
1;
