#ifndef SIMPLE_PSW_INCLUDED
#define SIMPLE_PSW_INCLUDED

#include "hmm/pairphmm.h"

/** Function to add references and acknowledgements text to an Opts_list class
 */
void add_simple_PSW_references (Opts_list& opts);

/** Parameters for simple Probabilistic Smith-Waterman */
struct Simple_PSW_params
{
  // alphabet
  const Alphabet& alphabet;
  // null model
  Prob null_extend;
  vector<Prob> null_emit;
  // gap open param for flanking regions
  Prob flank_open;
  // gap open/extend params within main alignment
  Prob gap_open;
  Prob gap_extend;
  // probability of another match, given that gap is finished
  Prob align_extend;
  // substitution matrix: submat(a,b) = probability of substituting b for a, conditional on having observed a
  array2d<Prob> submat;
  // constructor
  Simple_PSW_params (const Alphabet& alphabet);
  // helpers
  array2d<Prob> joint_submat() const;  // create joint substitution matrix from conditional matrix & null model
  // IO
  void write_submat (ostream& out) const;
  void read_submat (istream& in);
  void write_params (ostream& out) const;  // writes all params, including submat
  void read_params (istream& in);  // reads all params, including submat
  // file format help method
  static bool format_help (Opts_list* ol);
};

/** PPSW_params uses PAM250 as default, and has built-in defaults for gap params */
struct PPSW_params : Simple_PSW_params
{
  PPSW_params();
};

/** DPSW_params uses Hasegawa (i/v=10,p(a)=p(c)=p(g)=p(t)=1/4,t=1) as default, and has built-in defaults for gap params */
struct DPSW_params : Simple_PSW_params
{
  DPSW_params();
};

/** Simple_PSW_pscores offers a symbolic representation of the parameters for PSW
    that can be used for automated Baum-Welch training,
    using DART's "PGroup" class (and "Boolean_group" which is a subclass of PGroup),
    as well as the values of the PSW parameters in DART Score form (1/1000's of a bit)
*/
struct Simple_PSW_pscores : PScores
{
  const Alphabet& alphabet;
  // PGroup representations of params
  Boolean_group null_extend, flank_open, gap_open, gap_extend, align_extend;
  Alphabet_group null_emit;
  vector<Alphabet_group> submat;
  // constructor
  Simple_PSW_pscores (Simple_PSW_params& params);
  // methods for passing param values back & forth between Simple_PSW_params object
  void init_from_params (const Simple_PSW_params& params);
  void update_params (Simple_PSW_params& params) const;
  // method for optimising from a set of Baum-Welch counts (does not touch null params)
  void optimise_from_counts (const PCounts& pcounts);
};

/** Simple PSW Pair HMM, using PFunc class to represent transition/emission funcs */
struct Simple_PSW_PHMM : Pair_PHMM
{
  // state indices; NB match state is last (this is purely for convenience; the ins/del states are all lumped together)
  enum { LeftIns=0, LeftDel=1, RightIns=2, RightDel=3, Ins=4, Del=5, Mat=6, TotalStates=7 };
  // constructor
  Simple_PSW_PHMM (const Simple_PSW_pscores& pscores, bool odds_ratio = TRUE);
};

#endif /* SIMPLE_PSW_INCLUDED */
