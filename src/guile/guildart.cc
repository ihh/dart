#include <stdio.h>
#include <libguile.h>

#include "guile/stockholm-type.h"
#include "guile/newick-type.h"
#include "ecfg/guile-ecfg.h"
#include "ontology/onto_sexpr.h"

static void inner_main (void *closure, int argc, char **argv)
{
  try
    {
      init_stockholm_type();
      init_newick_type();

      init_dart_primitives();
      init_xrate_primitives();
      init_termx_primitives();

      scm_shell (argc, argv);
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }
}

int main (int argc, char **argv)
{
  SExpr_Scheme_evaluator::mark_guile_initialized();  // hack to avoid initializing guile twice
  scm_boot_guile (argc, argv, inner_main, 0);
  printf ("Guile exited - try installing guile and rebuilding\n");
  return 0;
}
