#include "psw/simplepsw.h"
#include "tree/pam.h"
#include "tree/hasegawa.h"
#include "util/vector_output.h"
#include "seq/dirichlet.h"

void add_simple_PSW_references (Opts_list& opts)
{
  opts.newline();
  opts.newline();
  opts.print ("Acknowledgements\n");
  opts.print ("----------------\n");
  opts.print ("Released under the GNU Public License by Ian Holmes, Dept of Bioengineering, UC Berkeley. Email <ihh@berkeley.edu>\n");
  opts.print ("Acknowledgements to users including: Michael Eisen, Colin Wiel, Robin Dowell.\n");
  opts.newline();
  opts.print ("References\n");
  opts.print ("----------\n");
  opts.print ("S.Miyazawa, 1994. A Reliable Sequence Alignment Method Based on Probabilities of Residue Correspondence. Protein Evolution, (8) 999-1009. [introduces the partition function for pairwise alignment]\n");
  opts.print ("P.Bucher and K.Hofmann, 1996. A Sequence Similarity Search Algorithm Based on a Probabilistic Interpretation of an Alignment Scoring System. Proc. 4th Intl Conf on Intelligent Systems for Molecular Biology, 44-51. [introduces Probabilistic Smith-Waterman using a partition function]\n");
  opts.print ("I.Holmes and R.Durbin, 1998. Dynamic Programming Alignment Accuracy. Journal of Computational Biology, 5(3) 493-504. [recasts sequence alignment as a Pair HMM; evaluates the optimal accuracy algorithm]\n");
  opts.print ("R.Durbin, S.Eddy, A.Krogh and G.Mitchison, 1998. Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press, Cambridge, UK. [general reference for Pair HMMs]\n");
}

Simple_PSW_params::Simple_PSW_params (const Alphabet& alphabet)
  : alphabet (alphabet),
    null_emit (alphabet.size(), 0.),
    flank_open(0.), gap_open(0.), gap_extend(0.), align_extend(0.),
    submat (alphabet.size(), alphabet.size(), 0.)
{ }

array2d<Prob> Simple_PSW_params::joint_submat() const
{
  array2d<Prob> joint = submat;
  for (int i = 0; i < alphabet.size(); ++i)
    for (int j = 0; j < alphabet.size(); ++j)
      joint(i,j) *= null_emit[i];
  return joint;
}

void Simple_PSW_params::write_submat (ostream& out) const
{
  out << null_emit << "\n";
  out << submat;
  submat.write_rowcol (out);
}

void Simple_PSW_params::read_submat (istream& in)
{
  for (int c = 0; c < alphabet.size(); ++c)
    in >> null_emit[c];
  submat.read_rowcol (in);
}

void Simple_PSW_params::write_params (ostream& out) const
{
  out << null_extend << " ";
  out << flank_open << " ";
  out << gap_open << " ";
  out << gap_extend << " ";
  out << align_extend << "\n";
  write_submat (out);
}

void Simple_PSW_params::read_params (istream& in)
{
  in >> null_extend;
  in >> flank_open;
  in >> gap_open;
  in >> gap_extend;
  in >> align_extend;
  read_submat (in);
}

PPSW_params::PPSW_params() : Simple_PSW_params (Protein_alphabet)
{
  // submat default: PAM250
  PAM_factory pam_factory;
  null_emit = pam_factory.create_prior();
  submat = pam_factory.create_conditional_substitution_matrix (2.50);
  // gap defaults
  null_extend = .6;
  flank_open = .5;
  gap_open = .01;
  gap_extend = .6;
  align_extend = .99;
}

DPSW_params::DPSW_params() : Simple_PSW_params (DNA_alphabet)
{
  // submat default: Hasegawa
  const double has_transition_rate = 1.0;
  const double has_transversion_rate = 1.0;
  const double has_time = .2;
  null_emit = vector<Prob> (4, 0.25);
  submat = Hasegawa_likelihood_matrix (has_transition_rate, has_transversion_rate, has_time, null_emit, FALSE, TRUE);
  // gap defaults
  null_extend = .6;
  flank_open = .5;
  gap_open = .01;
  gap_extend = .6;
  align_extend = .99;
}

Simple_PSW_pscores::Simple_PSW_pscores (Simple_PSW_params& params)
  : PScores(),
    alphabet (params.alphabet),
    null_extend (new_boolean_group ("nullExtend")),
    flank_open (new_boolean_group ("flankOpen")),
    gap_open (new_boolean_group ("gapOpen")),
    gap_extend (new_boolean_group ("gapExtend")),
    align_extend (new_boolean_group ("alignExtend")),
    null_emit (new_alphabet_group (params.alphabet, "nullEmit")),
    submat (params.alphabet.size())
{
  for (int c = 0; c < params.alphabet.size(); ++c)
    {
      sstring name;
      name << "subst" << params.alphabet.int2char_uc(c);
      submat[c] = new_alphabet_group (params.alphabet, name.c_str());
    }
  init_from_params (params);
}

void Simple_PSW_pscores::init_from_params (const Simple_PSW_params& params)
{
  PScores& sc (*this);

  sc[null_extend.y] = Prob2Score (params.null_extend);
  sc[null_extend.n] = Prob2Score (1. - params.null_extend);

  sc[flank_open.y] = Prob2Score (params.flank_open);
  sc[flank_open.n] = Prob2Score (1. - params.flank_open);

  sc[gap_open.y] = Prob2Score (params.gap_open);
  sc[gap_open.n] = Prob2Score (1. - params.gap_open);

  sc[gap_extend.y] = Prob2Score (params.gap_extend);
  sc[gap_extend.n] = Prob2Score (1. - params.gap_extend);

  sc[align_extend.y] = Prob2Score (params.align_extend);
  sc[align_extend.n] = Prob2Score (1. - params.align_extend);

  for (int c = 0; c < params.alphabet.size(); ++c)
    {
      sc[null_emit[c]] = Prob2Score (params.null_emit[c]);
      for (int d = 0; d < params.alphabet.size(); ++d)
	sc[submat[c][d]] = Prob2Score (params.submat(c,d));
    }
}

void Simple_PSW_pscores::update_params (Simple_PSW_params& params) const
{
  const PScores& sc (*this);

  params.null_extend = FScore2Prob (sc[null_extend.y]);
  params.flank_open = FScore2Prob (sc[flank_open.y]);
  params.gap_open = FScore2Prob (sc[gap_open.y]);
  params.gap_extend = FScore2Prob (sc[gap_extend.y]);
  params.align_extend = FScore2Prob (sc[align_extend.y]);

  for (int c = 0; c < params.alphabet.size(); ++c)
    {
      params.null_emit[c] = FScore2Prob (sc[null_emit[c]]);
      for (int d = 0; d < params.alphabet.size(); ++d)
	params.submat(c,d) = FScore2Prob (sc[submat[c][d]]);
    }
}

void Simple_PSW_pscores::optimise_from_counts (const PCounts& pcounts)
{
  Multigroup multi = flank_open * gap_open * gap_extend * align_extend * null_extend * null_emit;
  for (int c = 0; c < alphabet.size(); ++c)
    multi *= submat[c];
  Laplace_prior laplace (multi);
  Dirichlet_assignment dir (laplace, multi);
  dir.optimise_scores (pcounts, *this);
}

Simple_PSW_PHMM::Simple_PSW_PHMM (const Simple_PSW_pscores& pscores, bool odds_ratio)
  : Pair_PHMM (TotalStates, pscores.alphabet)
{
  // make references to PGroups
  const Boolean_group& null_extend = pscores.null_extend;
  const Boolean_group& flank_open = pscores.flank_open;
  const Boolean_group& flank_extend = pscores.null_extend;  // NB flank_extend == null_extend
  const Boolean_group& gap_open = pscores.gap_open;
  const Boolean_group& gap_extend = pscores.gap_extend;
  const Boolean_group& align_extend = pscores.align_extend;
  const Alphabet_group& null_emit = pscores.null_emit;
  const vector<Alphabet_group>& submat = pscores.submat;
  // set state types
  state_type[LeftIns] = state_type[Ins] = state_type[RightIns] = EmitY;
  state_type[LeftDel] = state_type[Del] = state_type[RightDel] = EmitX;
  state_type[Mat] = EmitXY;
  // set single emission scores
  for (int s = LeftIns; s < Mat; ++s)  // loop through all ins/del states
    {
      single_emit[s] = vector<PFunc> (alphabet->size());
      for (int c = 0; c < alphabet->size(); ++c)
	if (odds_ratio)
	  single_emit[s][c] = 1. / null_extend.y;
	else
	  single_emit[s][c] = null_emit[c];
    }
  // set pair emission scores
  pair_emit[Mat] = array2d<PFunc> (alphabet->size(), alphabet->size());
  for (int c = 0; c < alphabet->size(); ++c)
    for (int d = 0; d < alphabet->size(); ++d)
      if (odds_ratio)
	pair_emit[Mat](c,d) = submat[c][d] / (null_emit[d] * null_extend.y * null_extend.y);
      else
	pair_emit[Mat](c,d) = null_emit[c] * submat[c][d];
  // set transitions
  // transitions from Start
  transition(Start,    LeftIns)  = flank_open.y;
  transition(Start,    LeftDel)  = flank_open.n * flank_open.y;
  transition(Start,    Mat)      = flank_open.n * flank_open.n * align_extend.y;
  transition(Start,    End)      = flank_open.n * flank_open.n * align_extend.n;
  // transitions from LeftIns
  transition(LeftIns,  LeftIns)  = flank_extend.y;
  transition(LeftIns,  LeftDel)  = flank_extend.n * flank_open.y;
  transition(LeftIns,  Mat)      = flank_extend.n * flank_open.n * align_extend.y;
  transition(LeftIns,  End)      = flank_extend.n * flank_open.n * align_extend.n;
  // transitions from LeftDel
  transition(LeftDel,  LeftDel)  = flank_extend.y;
  transition(LeftDel,  Mat)      = flank_extend.n * align_extend.y;
  transition(LeftDel,  End)      = flank_extend.n * align_extend.n;
  // transitions from Mat
  transition(Mat,      Ins)      = align_extend.y * gap_open.y;
  transition(Mat,      Del)      = align_extend.y * gap_open.n * gap_open.y;
  transition(Mat,      Mat)      = align_extend.y * gap_open.n * gap_open.n;
  transition(Mat,      RightIns) = align_extend.n * flank_open.y;
  transition(Mat,      RightDel) = align_extend.n * flank_open.n * flank_open.y;
  transition(Mat,      End)      = align_extend.n * flank_open.n * flank_open.n;
  // transitions from Ins
  transition(Ins,      Ins)      = gap_extend.y;
  transition(Ins,      Del)      = gap_extend.n * gap_open.y;
  transition(Ins,      Mat)      = gap_extend.n * gap_open.n;
  // transitions from Del
  transition(Del,      Del)      = gap_extend.y;
  transition(Del,      Mat)      = gap_extend.n;
  // transitions from RightIns
  transition(RightIns, RightIns) = flank_extend.y;
  transition(RightIns, RightDel) = flank_extend.n * flank_open.y;
  transition(RightIns, End)      = flank_extend.n * flank_open.n;
  // transitions from RightDel
  transition(RightDel, RightDel) = flank_extend.y;
  transition(RightDel, End)      = flank_extend.n;
  // if using odds ratios, divide end transitions by (null_extend.n)^2
  if (odds_ratio)
    for (int src = Start; src < TotalStates; ++src)
      if (src != End)
	transition(src,End) /= null_extend.n * null_extend.n;
}

bool Simple_PSW_params::format_help (Opts_list* ol)
{
  CLOGERR << "\n";
  CLOGERR << "File format for PSW substitution matrix\n";
  CLOGERR << "=======================================\n";
  CLOGERR << "Suppose A is the alphabet size; then there are A+1 lines in the file.\n";
  CLOGERR << "The first line has A columns; column i is P(i), the initial probability of residue i.\n";
  CLOGERR << "The subsequent (A) lines each have (A) fields; row i column j is P(j|i), the probability of j conditional on i.\n";
  CLOGERR << "All fields are delimited by whitespace.\n";
  CLOGERR << "\n";
  CLOGERR << "File format for extended PSW parameters\n";
  CLOGERR << "=======================================\n";
  CLOGERR << "The first line contains the parameters (null_extend, flank_open, gap_open, gap_extend, align_extend) in that order.\n";
  CLOGERR << "The subsequent lines contain the PSW substitution matrix (see above).\n";
  exit(0);
  return 0;
}
