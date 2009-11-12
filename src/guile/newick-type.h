#ifndef NEWICK_TYPE_INCLUDED
#define NEWICK_TYPE_INCLUDED

#include "guile/guile-defs.h"
#include "tree/phylogeny.h"

// guile functions
extern scm_t_bits newick_tag;
void init_newick_type (void);
SCM make_newick_smob (const PHYLIP_tree& tree);

// guile smobs
// cast method
PHYLIP_tree* newick_cast_from_scm (SCM tree_smob);

#endif /* NEWICK_TYPE_INCLUDED */
