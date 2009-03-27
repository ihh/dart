#ifndef ECFG_MAIN_INCLUDED
#define ECFG_MAIN_INCLUDED

#include <sstream>
#include "util/logfile.h"
#include "util/rnd.h"
#include "util/vector_output.h"
#include "tree/subdistmat.h"
#include "ecfg/pfold.h"
#include "ecfg/ecfgdp.h"

struct ECFG_main
{
  // typedefs
  typedef map<sstring,ECFG_scores*> ECFG_map;

  // data
  Opts_list opts;
  const Alphabet* alph;
  ECFG_chain* tree_estimation_chain;  // contains matrix for doing neighbor-joining, branch-length EM
  ECFG_map ecfg_map;  // list of grammars

  // command-line params
  bool annotate;
  bool ancrec_CYK_MAP;
  bool ancrec_postprob;
  bool report_sumscore;
  bool report_confidence;
  bool report_postprob;
  bool report_hidden_classes;
  int max_subseq_len;
  sstring grammars_filename;
  sstring tree_grammar_filename;
  sstring dump_expanded;
  sstring preset;
  sstring train;
  sstring gff_filename;
  sstring gap_chars;

  double em_min_inc;
  int em_max_iter;
  int em_forgive;
  double tres;

  double pseud_init;
  double pseud_mutate;
  double pseud_wait;

  double min_branch_len;

  vector<sstring> grammar_list;
  sstring default_grammars;

  vector<ECFG_scores*> grammar;  // actual ECFG objects

  // constructor
  ECFG_main (int argc, char** argv);

  // method to add grammar
  void add_grammar (const char* name, ECFG_scores* ecfg);

  // method to annotate a score
  void annotate_loglike (Stockholm& stock, const char* tag, const sstring& ecfg_name, Loge loglike) const;

  // opts-list initialiser
  void init_opts (const char* desc);

  // method
  void run();
};

#endif /* ECFG_MAIN_INCLUDED */
