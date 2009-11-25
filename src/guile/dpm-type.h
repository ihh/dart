#ifndef XRATE_DPM_TYPE_INCLUDED
#define XRATE_DPM_TYPE_INCLUDED

#include "guile/guile-defs.h"
#include "ecfg/ecfgdp.h"

// guile functions
extern scm_t_bits xrate_dpm_tag;
void init_xrate_dpm_type (void);
SCM make_xrate_dpm_smob (const ECFG_scores& ecfg, const Stockholm& stock);

// guile smobs
struct XRATE_DPM_smob {
  // Scheme objects
  SCM stock_smob_scm;
  // internal data
  ECFG_scores* ecfg;
  ECFG_envelope* env;
  PCounts* pcounts;
  ECFG_inside_outside_matrix* inout;
  ECFG_CYK_matrix* cyk;

  // constructors
  XRATE_DPM_smob (SCM alphabet_and_grammar_scm, SCM stock_smob_scm);

  // destructor
  ~XRATE_DPM_smob();

  // cast method
  static XRATE_DPM_smob* cast_from_scm (SCM xrate_dpm_smob);
};

#endif /* XRATE_DPM_TYPE_INCLUDED */
