#include <libguile.h>
#include "guile/stockholm-type.h"

static void inner_main (void *closure, int argc, char **argv)
{
  init_stockholm_type();
  init_xrate_primitives();
  scm_shell (argc, argv);
}

int main (int argc, char **argv)
{
  scm_boot_guile (argc, argv, inner_main, 0);
  return 0; /* never reached */
}
