#ifndef ECFG_MAIN_INCLUDED
#define ECFG_MAIN_INCLUDED

#include <sstream>
#include "util/logfile.h"
#include "util/rnd.h"
#include "util/vector_output.h"
#include "tree/subdistmat.h"
#include "hsm/branch_length_em.h"
#include "ecfg/pfold.h"
#include "ecfg/ecfgdp.h"

struct ECFG_main
{
  // typedefs
  typedef map<sstring,ECFG_scores*> ECFG_map;

  // data
  // command-line option parser
  Opts_list opts;

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

  // alphabet & grammar data
  const Alphabet* alph;
  Empty_alphabet user_alphabet;

  sstring default_grammars;
  vector<sstring> grammar_list; // grammar names
  vector<ECFG_scores*> grammar;  // grammar objects
  ECFG_map ecfg_map;  // map from grammar names to grammar objects

  const ECFG_chain* tree_estimation_chain;  // contains matrix for doing neighbor-joining, branch-length EM

  // alignment data
  Sequence_database seq_db;
  Stockholm_database stock_db;
  vector<sstring> training_alignment_filename;
  EM_tree_alignment_database align_db;
  vector<Aligned_score_profile> asp_vec;

  // constructors
  ECFG_main (int argc, char** argv);

  // add a preset grammar
  // call this as many times as you like after constructing the object
  void add_grammar (const char* name, ECFG_scores* ecfg);

  // initialise opts-list
  // call this after constructing & adding preset grammars
  void init_opts (const char* desc);

  // top-level run method for xrate
  void run (ostream& alignment_output_stream);

  // alternate top-level methods for embedded invocation (currently placeholders)
  Stockholm estimate_tree(const ECFG_scores& ecfg, Stockholm& stock);
  ECFG_scores train_grammar(const ECFG_scores& ecfg, Stockholm_database& stock);
  Stockholm annotate_alignment(const ECFG_scores& ecfg, Stockholm& stock);

  // lower-level methods called by run()
  // comments indicate required modifications for embedded invocation
  void parse_opts();  // can be skipped
  void read_alignments();  // must be called; add optional argument to allow direct specification of an alignment?
  void estimate_trees();  // can be skipped; add optional argument to allow direct specification of a tree-estimation grammar?
  void read_grammars();  // must be called (but can use a default grammar); add optional argument to allow direct specification of a grammar?
  void convert_sequences();  // must be called if train_grammars() or annotate_alignments() is to be called
  void train_grammars();  // can be skipped
  void annotate_alignments (ostream* align_stream = 0);  // can be skipped

  // helper to add a particular score annotation to an alignment
  void annotate_loglike (Stockholm& stock, const char* tag, const sstring& ecfg_name, Loge loglike) const;

  // helper to test if any trees are missing
  bool missing_trees() const;
};

#endif /* ECFG_MAIN_INCLUDED */
