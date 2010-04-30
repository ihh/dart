#include <stdio.h>

#ifdef GUILE_INCLUDED
#include <libguile.h>
#include "guile/stockholm-type.h"
#include "guile/newick-type.h"
#include "ecfg/guile-ecfg.h"
#endif /* GUILE_INCLUDED */

#ifdef GUILE_INCLUDED
static void inner_main (void *closure, int argc, char **argv)
{
  init_stockholm_type();
  init_newick_type();
  init_xrate_primitives();
  scm_shell (argc, argv);
}
#endif /* GUILE_INCLUDED */

int main (int argc, char **argv)
{
#ifdef GUILE_INCLUDED
  scm_boot_guile (argc, argv, inner_main, 0);
#endif /* GUILE_INCLUDED */
  printf ("Guile unavailable - try installing guile and rebuilding\n");
  return 0; /* only reached if GUILE_INCLUDED is not defined */
}
