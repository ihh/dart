# Perl module for working with xgram grammar files in protein sequences

package PhyloGram::Protein;

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

    # add Protein alphabet
    $self->add (DartSexpr->from_string (<<END));
alphabet
(name Protein)
(token (a r n d c q e g h i l k m f p s t w y v))
(extend (to x) (from a) (from r) (from n) (from d) (from c) (from q) (from e) (from g) (from h) (from i) (from l) (from k) (from m) (from f) (from p) (from s) (from t) (from w) (from y) (from v))
(extend (to b) (from n) (from d))
(extend (to z) (from q) (from e))
(wildcard *)
END

    # return
    return $self;
}

# end of package
1;
