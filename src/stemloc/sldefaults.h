#ifndef SLDEFAULTS_INCLUDED
#define SLDEFAULTS_INCLUDED

#include "telegraph/tgvar.h"
#include "scfg/cfgsexpr.h"

struct Stemloc_defaults
{
  static void init_quickstem_params (Telegraph_PScores_adaptor& fold_tgio);
  static void init_quickalign_params (Telegraph_PScores_adaptor& align_tgio);
  static void init_stempair_params (Telegraph_PScores_adaptor& cmp_tgio);

  static void init_superstem_params (Telegraph_PScores_adaptor& fold_tgio);
  static void init_superpair_params (Telegraph_PScores_adaptor& cmp_tgio);

  static void init_flat_params (Telegraph_PScores_adaptor& tgio);
};

class Stemloc_gramset : public Gramset
{
private:
  // helper subroutine to "fix" group_suffix for a particular PScores
  void fix_gs (PScores& pscores);
public:
  Stemloc_gramset (bool fix_group_suffix_now = true);  // there is really no need to leave fix_group_suffix_now unset, except to convert old Telegraph-format files to new stemloc format
  void fix_group_suffix();  // called by constructor if appropriate flag is set
  static void prefix_comments (ostream& out);  // dumps some explanatory comments on model format to ostream...  in lieu of proper documentation... :-&
};

#endif /* SLDEFAULTS_INCLUDED */
