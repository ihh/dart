// pk 5/05 test xfold using Pair_CFG with pfold's grammar, parameters - based on evoldoer.cc, testpaircfgdp.cc
// 8/05 use PFOLD_loop_matrix, PFOLD_stem_matrix

#include <sstream>
#include "scfg/postenv.h"
#include "scfg/cfgdotplot.h"
#include "paircfg_xfold.h"

// main
int main (int argc, char** argv)
{
  try
  {
    // create default parameters
    PFOLD_TKFST_params params;

    // create opts
    INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <alignment database in Stockholm format>",
		"generate inside score for alignment using Pair_CFG\n");

    opts.newline();
    double t;
    double half;
    bool CYK;
    opts.add ("t -time", t = 1, "evolutionary time");
    opts.add ("ht -ht half -half", half = -1, "half of the evolutionary time; overrides --time", false);
    opts.add ("cyk -cyk CYK -CYK", CYK = 1, "use CYK matrix only");

	sstring pfold_prefix;
	opts.add ("pff -postfoldfile", pfold_prefix = "", "filename prefix for fold post.prob. matrices", FALSE);

    // parse opts
    opts.parse_or_die();
    if (half >= 0)	// half overrides time
	  t = half * 2;

    // get alignment filename
    const char* alignment_db_filename = opts.args[0].c_str();

    // create databases & alphabet
    const Alphabet& alphabet (CFG_alphabet);
    Stockholm_database align_db;
    FASTA_sequence_database seq_db;

    ifstream align_db_in (alignment_db_filename);
    if (!align_db_in) THROWEXPR ("Couldn't open alignment file '" << alignment_db_filename << "'");
    align_db.read_Stockholm_or_FASTA (align_db_in, seq_db);

    // update the sequence database index
    seq_db.update_index();
    seq_db.seqs_update (alphabet, (Profile_flags_enum::Profile_flags) (Profile_flags_enum::DSQ | Profile_flags_enum::SCORE));
      
    // create a PFOLD Pair SCFG
    const PFOLD_pair_CFG cfg (params, t);

    // need at least two sequences - use seq 0 for x, seq 1 for y
    if (seq_db.size() < 2)
	  THROWEXPR("Need at least two sequences for test");

    // follow test_pair_emit()
    const Named_profile& npx = *seq_db.index.profile[0];
    const Named_profile& npy = *seq_db.index.profile[1];

    Fold_envelope xenv;
    Fold_envelope yenv;

	if (npx.dsq.size() != npy.dsq.size())	// sequences must have the same lenght (non-gap)
	  THROWEXPR ("Sequence mismatch: x size " << npx.dsq.size() << " y size " << npy.dsq.size());
    xenv.initialise_full (npx.dsq.size());
    yenv.initialise_full (npy.dsq.size());

    // set up Pairwise_path
    const Stockholm_database::Align_row& arx = align_db.align_row[npx.name];
    const Stockholm_database::Align_row& ary = align_db.align_row[npy.name];
    if (arx.nAlign != ary.nAlign)
	  THROWEXPR("Invalid alignment");
    const Pairwise_path path (align_db.align_index[arx.nAlign]->path, arx.nRow, ary.nRow, TRUE);
    Pair_envelope pair_env (npx.size(), npy.size(), 0);
    pair_env.add_pairwise_path (path);

    // calculate, print desired matrix
    if (CYK)
    {
	  Pair_CYK_matrix cyk (npx, npy, xenv, yenv, cfg, pair_env);
	  cout << "Pair_CYK_matrix score " << Score2Bits (cyk.final_score) << "\n";
	  // get alignment
	  const Pair_CFG_alignment cyk_alignment = cyk.alignment();
	  // output
	  cyk_alignment.show (cout);
    }
	else if (pfold_prefix.size())
	{
	  const Pair_inside_outside_matrix in_out (npx, npy, xenv, yenv, cfg, pair_env, FALSE);
	  const PairCFG_fold_dotplot xfold_dotplot (in_out, 0);
	  xfold_dotplot.write_dotplot (pfold_prefix, npx.name);

	  const PairCFG_fold_dotplot yfold_dotplot (in_out, 1);
	  yfold_dotplot.write_dotplot (pfold_prefix, npy.name);
	}
    else
    {
	  Pair_inside_matrix inside (npx, npy, xenv, yenv, cfg, pair_env);
	  cout << "Pair_inside_matrix score " << Score2Bits (inside.final_score) << "\n";
    }

  }		// end try
  catch (const Dart_exception& e)
  {
	CLOGERR << e.what();
	return 1;
  }
  return 0;
}

#define INIT_LAMBDA 1.
#define INIT_MU     2.

PFOLD_TKFST_params::PFOLD_TKFST_params()
  : loop (loop_mat, INIT_LAMBDA, INIT_MU),
    stem (stem_mat, INIT_LAMBDA, INIT_MU)
{
  // Stem & loop indel params (almost totally arbitrary)
  stem.lambda = 0.007;
  stem.mu     = 0.01;
  loop.lambda = 0.027;
  loop.mu     = 0.03;
  stem_prob   = 0.01;
}

PFOLD_pair_CFG::PFOLD_pair_CFG (PFOLD_TKFST_params& params, double time)
  : Pair_CFG_scores (pfoldTotalStates),
    params (params),
    time (time)
{
  // set state types
  init_emit (pfoldU, EmitXLYL);
  init_emit (pfoldF, EmitXLRYLR);

  state_type[pfoldS]       = Null;
  state_type[pfoldL]       = Null;

  state_type[pfoldB]      = Bifurc;

  // make bifurcations
  bifurc[pfoldB]      = Bifurcation (pfoldL, pfoldS);

  // calculate some score parameters
  const array2d<Score> loop_pair_sc = Prob2ScoreArray2d (params.loop_mat.create_joint_substitution_matrix (time));
  const array2d<Score> stem_pair_sc = Prob2ScoreArray2d (params.stem_mat.create_joint_substitution_matrix (time));

  // make emissions
  for (int xl = 0; xl < CFG_alphabet_size; ++xl)
  {
	for (int yl = 0; yl < CFG_alphabet_size; ++yl)
	{
	  const int xlyl_idx = emit_xl_mul(EmitXLYL) * xl + emit_yl_mul(EmitXLYL) * yl;
	  emit[pfoldU][xlyl_idx] = loop_pair_sc(xl,yl);
	}
	for (int xr = 0; xr < CFG_alphabet_size; ++xr)
	{
	  const int xlr_idx = emit_xl_mul(EmitXLR) * xl + emit_xr_mul(EmitXLR) * xr;
	  for (int yl = 0; yl < CFG_alphabet_size; ++yl)
	    for (int yr = 0; yr < CFG_alphabet_size; ++yr)
		{
		  const int ylr_idx = emit_yl_mul(EmitYLR) * yl + emit_yr_mul(EmitYLR) * yr;
		  const int xlrylr_idx = emit_xl_mul(EmitXLRYLR) * xl + emit_xr_mul(EmitXLRYLR) * xr + emit_yl_mul(EmitXLRYLR) * yl + emit_yr_mul(EmitXLRYLR) * yr;
		  emit[pfoldF][xlrylr_idx] = stem_pair_sc(xlr_idx,ylr_idx);
		}
	}
  }

  // make transitions
  transition (Start, pfoldS) = Prob2Score (1);
  transition (pfoldS, pfoldB) = Prob2Score (0.868534);
  transition (pfoldS, pfoldL) = Prob2Score (0.131466);
  transition (pfoldL, pfoldU) = Prob2Score (0.894603);
  transition (pfoldL, pfoldF) = Prob2Score (0.105397);
  transition (pfoldU, End) = Prob2Score (1);
  transition (pfoldF, pfoldF) = Prob2Score (0.787640);
  transition (pfoldF, pfoldB) = Prob2Score (0.212360);

}
