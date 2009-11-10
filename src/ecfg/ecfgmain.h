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

  // training data
  vector<ECFG_trainer*> trainer;

  // constructors
  ECFG_main (int argc, char** argv);  // does not call init_opts, so caller can add grammars before initializing the Opts_List
  ECFG_main();  // calls init_opts implicitly, for embedded usage

  // destructor
  ~ECFG_main() { delete_trainers(); }

  // add a preset grammar
  // call this as many times as you like after constructing the object
  void add_grammar (const char* name, ECFG_scores* ecfg);

  // initialise opts-list
  // call this after constructing & adding preset grammars
  void init_opts (const char* desc);

  // top-level run method for xrate
  // note that none of the run methods catch exceptions...
  void run_xrate (ostream& alignment_output_stream);

  // alternate top-level run methods for embedded invocation (currently placeholders)
  // note that none of the run methods catch exceptions...
  Stockholm& run_tree_estimation (Stockholm& stock, Sequence_database& seq_db, SExpr& grammar_alphabet_sexpr);
  Stockholm& run_alignment_annotation (Stockholm& stock, SExpr& grammar_alphabet_sexpr);
  void run_grammar_training (Stockholm_database& stock, SExpr& grammar_alphabet_sexpr, ECFG_scores** grammar_ret, ECFG_counts** counts_ret);

  // lower-level methods called by run_*()
  // comments indicate required usage for embedded invocation
  void parse_opts();  // can be skipped
  void read_alignments (const Stockholm_database* stock = 0);  // must be called; optional argument allows direct specification of an alignment
  void read_alignments (const Stockholm& stock);  // overloaded form of read_alignments allowing specification of a single alignment
  void estimate_trees (SExpr* grammar_alphabet_sexpr = 0, Sequence_database* seq_db_ptr = 0);  // can be skipped; optional arguments allow direct specification of a tree-estimation grammar and sequence database
  void read_grammars (SExpr* grammar_alphabet_sexpr = 0);  // must be called (but can use a default grammar); optional argument allows direct specification of an alphabet/grammar pair
  void convert_sequences();  // must be called if train_grammars() or annotate_alignments() is to be called
  void train_grammars();  // can be skipped
  void delete_trainers();  // can be skipped if train_grammars() was skipped
  void annotate_alignments (ostream* align_stream = 0);  // can be skipped

  // helper to add a particular score annotation to an alignment
  void annotate_loglike (Stockholm& stock, const char* tag, const sstring& ecfg_name, Loge loglike) const;

  // helper to test if any trees are missing
  bool missing_trees() const;
};

#endif /* ECFG_MAIN_INCLUDED */
