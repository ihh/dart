#include "tkf/tkfopts.h"
#include "tree/hasegawa.h"
#include "tree/pam.h"
#include "util/logfile.h"
#include "hsm/em_matrix.h"
#include "util/vector_output.h"

TKF_opts::TKF_opts (int argc, char** argv, bool input) :
  Opts_list (argc, argv),
  bd_matrix(2,2),
  bd_vec(2),
  iv_matrix(2,2),
  iv_vec(2),
  composition(4,.25),
  submat_factory(0)
{
  SET_VERSION_INFO (*this);
  Rnd::add_opts (*this);
  newline();
  Log_stream::add_opts (*this);

  print ("Indel model parameters (Thorne et al, 1991)\n");
  print ("-------------------------------------------\n");

  bd_matrix(1,1) = 1.02; bd_matrix(1,2) = -1; bd_vec(1) = 0;   last_bd_constraint = SEQLEN;
  bd_matrix(2,1) = 1;    bd_matrix(2,2) = 1;  bd_vec(2) = 0.1; last_bd_constraint = INDEL;

  iv_matrix(1,1) = 1;    iv_matrix(1,2) = 2;  iv_vec(1) = 1;   last_iv_constraint = SUBST;
  iv_matrix(2,1) = 1;    iv_matrix(2,2) = -4; iv_vec(2) = 0;   last_iv_constraint = IVRATIO;

  add ("seqlen",    &specify_bd_seqlen,    "\texpected sequence length }");
  add ("indelrate", &specify_bd_indelrate, "\tindel rate               } need two of these to specify");
  add ("birthrate", &specify_bd_birthrate, "\tbirth rate               } birth and death rates");
  add ("deathrate", &specify_bd_deathrate, "\tdeath rate               } (default is -seqlen 50 -indelrate 0.1)");

  newline();
  print ("DNA substitution model parameters (Hasegawa et al, JME 1985)\n");
  print ("------------------------------------------------------------\n");

  add ("irate",   &specify_iv_irate,   "\ttransition rate               }");
  add ("vrate",   &specify_iv_vrate,   "\ttransversion rate             } need two of these to specify");
  add ("ivratio", &specify_iv_ivratio, "\ttransition/transversion ratio } transition and transversion rates");
  add ("subrate", &specify_iv_subst,   "\tsubstitution rate             } (default is -subrate 1 -ivratio 4)");
  add ("gccontent",   &specify_gc_content,       " <float>\tspecify G+C content");
  add ("atcontent",   &specify_at_content,       " <float>\tspecify A+T content");
  add ("acgt",        &specify_composition,      " <a> <c> <g> <t>\tspecify full DNA composition vector");
  add ("composition", &specify_composition_file, " <string>\tread DNA composition vector from file");

  newline();
  print ("Alternative substitution models\n");
  print ("-------------------------------\n");
  
  add ("hsm", hsm_filename,        "read 'xrate'-format substitution model from file (-hsmhelp for format description)", 0);
  add ("pam", use_pam = 0, "use PAM substitution rate matrix");
  if (input)
    add ("autoseq", auto_alph = 1, "auto-detect nucleotide (Hasegawa) or protein (PAM)");
  add ("hsmhelp", &EM_matrix::hsm_help);

  post_parse_callback.push_back (&log_params);
}

bool TKF_opts::specify_bd_constraint (Opts_list* ol, BD_constraint constraint_type, double b_coeff, double d_coeff, double result)
{
  TKF_opts* t = (TKF_opts*) ol;
  if (constraint_type != t->last_bd_constraint)
    {
      t->bd_matrix.Row(1) = t->bd_matrix.Row(2);
      t->bd_vec(1) = t->bd_vec(2);
      t->last_bd_constraint = constraint_type;
    }
  
  t->bd_matrix(2,1) = b_coeff;
  t->bd_matrix(2,2) = d_coeff;
  t->bd_vec(2) = result;

  return 1;
}

bool TKF_opts::specify_iv_constraint (Opts_list* ol, IV_constraint constraint_type, double i_coeff, double v_coeff, double result)
{
  TKF_opts* t = (TKF_opts*) ol;
  if (constraint_type != t->last_iv_constraint)
    {
      t->iv_matrix.Row(1) = t->iv_matrix.Row(2);
      t->iv_vec(1) = t->iv_vec(2);
      t->last_iv_constraint = constraint_type;
    }
  
  t->iv_matrix(2,1) = i_coeff;
  t->iv_matrix(2,2) = v_coeff;
  t->iv_vec(2) = result;

  return 1;
}

bool TKF_opts::log_params (Opts_list* ol)
{
  TKF_opts* tkf = (TKF_opts*) ol;

  ColumnVector bd = tkf->get_bd();
  CLOG(8) << "Insertion rate = " << bd(1) << ", deletion rate = " << bd(2) << "\n";

  if (tkf->hsm_filename != "") CLOG(8) << "Hidden substitution matrix in file " << tkf->hsm_filename << "\n";
  else if (tkf->use_pam)
    {
      CLOG(8) << "Using PAM matrix\n";
    }
  else
    {
      ColumnVector iv = tkf->get_iv();
      CLOG(8) << "Transition rate = " << iv(1) << ", transversion rate = " << iv(2) << "\n";
      CLOG(8) << "Base frequency vector = (" << tkf->composition << ")\n";
    }

  return 1;
}

void TKF_opts::detect_alphabet (const Sequence_database& seq_db)
{
  if (auto_alph)
    use_pam = use_pam || (&seq_db.detect_alphabet() == &Protein_alphabet);
}

TKF_params TKF_opts::params()
{
  if (hsm_filename != "")
    {
      ifstream in (hsm_filename.c_str());
      if (!in) THROW Standard_exception ("Hidden subsitution matrix file not found");
      EM_matrix* em_matrix = new EM_matrix (2, 2, 1);  // dummy params
      em_matrix->read (in);
      em_matrix->guess_alphabet();
      submat_factory = em_matrix;
    }
  else if (use_pam)
    {
      submat_factory = new PAM_factory();
    }
  else
    {
      ColumnVector iv = get_iv();
      submat_factory = new Hasegawa_matrix_factory (iv(1), iv(2), composition);
    }

  ColumnVector bd = get_bd();
  if (bd(1) == 0 || bd(2) == 0)
    THROWEXPR ("Can't have zero insertion/deletion rate; use e.g. 0.0000001 instead");

  return TKF_params (*submat_factory, bd(1), bd(2));
}

void TKF_opts::add_align_detection()
{
  expect_aligned = AUTO;
  add ("aligned", &specify_aligned, " yes|no|auto\tuse pre-aligned input sequences (default is auto)");
}

bool TKF_opts::detect_aligned (istream& in) const
{
  if (expect_aligned == NO) return 0;
  if (expect_aligned == YES) return 1;
  return !Named_profile::detect_FASTA(in);
}

