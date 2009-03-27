#ifndef TKFOPTS_INCLUDED
#define TKFOPTS_INCLUDED

#include "tkf/tkfparams.h"
#include "tree/substitution_matrix_factory.h"
#include "util/opts_list.h"

struct TKF_opts : Opts_list
{
  // to specify birth & death rates requires two (out of four) independent linear constraints

  enum BD_constraint { SEQLEN, INDEL, BIRTH, DEATH };

  Matrix        bd_matrix;
  ColumnVector  bd_vec;
  BD_constraint last_bd_constraint;

  static bool specify_bd_constraint (Opts_list* ol, BD_constraint constraint_type, double b_coeff, double d_coeff, double result);

  static bool specify_bd_seqlen (Opts_list* ol)
    { return specify_bd_constraint (ol, SEQLEN, 1 + 1 / ol->next_double(), -1, 0); }

  static bool specify_bd_indelrate (Opts_list* ol)
    { return specify_bd_constraint (ol, INDEL, 1, 1, ol->next_double()); }

  static bool specify_bd_birthrate (Opts_list* ol)
    { return specify_bd_constraint (ol, BIRTH, 1, 0, ol->next_double()); }
 
  static bool specify_bd_deathrate (Opts_list* ol)
    { return specify_bd_constraint (ol, DEATH, 0, 1, ol->next_double()); }

  ColumnVector get_bd() { return bd_matrix.i() * bd_vec; }
  
  // rate matrix for hidden substitution model
  sstring hsm_filename;

  // PAM flag & auto-alphabet detect

  bool use_pam;
  bool auto_alph;

  // hasegawa params
  //
  // to specify hasegawa transition & transversion rates requires two (out of four) independent linear constraints

  enum IV_constraint { SUBST, IVRATIO, IRATE, VRATE };

  Matrix        iv_matrix;
  ColumnVector  iv_vec;
  IV_constraint last_iv_constraint;

  static bool specify_iv_constraint (Opts_list* ol, IV_constraint constraint_type, double i_coeff, double v_coeff, double result);

  static bool specify_iv_subst (Opts_list* ol)
    { return specify_iv_constraint (ol, SUBST, 1, 2, ol->next_double()); }

  static bool specify_iv_ivratio (Opts_list* ol)
    { return specify_iv_constraint (ol, IVRATIO, 1, -ol->next_double(), 0); }

  static bool specify_iv_irate (Opts_list* ol)
    { return specify_iv_constraint (ol, IRATE, 1, 0, ol->next_double()); }
 
  static bool specify_iv_vrate (Opts_list* ol)
    { return specify_iv_constraint (ol, VRATE, 0, 1, ol->next_double()); }

  ColumnVector get_iv() { return iv_matrix.i() * iv_vec; }

  // various ways of specifying residue frequencies for hasegawa model
  
  vector<double> composition;

  static bool specify_gc_content (Opts_list* ol)
    {
      TKF_opts* t = (TKF_opts*) ol;
      double gc = ol->next_double();
      t->composition[0] = t->composition[2] = gc/2;
      t->composition[1] = t->composition[3] = (1-gc)/2;
      return 1;
    }

  static bool specify_at_content (Opts_list* ol)
    {
      TKF_opts* t = (TKF_opts*) ol;
      double at = ol->next_double();
      t->composition[1] = t->composition[3] = at/2;
      t->composition[0] = t->composition[2] = (1-at)/2;
      return 1;
    }

  static bool specify_composition (Opts_list* ol)
    { for (int i = 0; i < 4; i++) ((TKF_opts*)ol)->composition[i] = ol->next_double(); return 1; }

  static bool specify_composition_file (Opts_list* ol)
    { ifstream in (ol->next_string()); for (int i = 0; i < 4; i++) in >> ((TKF_opts*)ol)->composition[i]; return 1; }

  Substitution_matrix_factory* submat_factory;

  Expect_flag expect_aligned;
  static bool specify_aligned (Opts_list* ol)
    { ol->set_expect_flag (((TKF_opts*)ol)->expect_aligned, ol->next_string()); return 1; }

  static bool log_params (Opts_list* ol);

  // member functions
  //
  TKF_opts (int argc, char** argv, bool input = 0);
  ~TKF_opts() { if (submat_factory) delete submat_factory; }

  void add_align_detection();
  bool detect_aligned (istream& in) const;

  void detect_alphabet (const Sequence_database& seq_db);  // sets use_pam flag

  TKF_params params();
};

#endif
