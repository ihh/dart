#include "scfg/cfgsexpr.h"
#include "stemloc/sldefaults.h"

// number of comparison parameterisations
#define N_CMP_PARAMS 4

// main program
int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, -1, "[options] [<stemloc gramset file>]",
		  "dump stemloc gramset to stdout");

  opts.print_title ("Legacy Telegraph format conversion options");

  sstring param_set, fold_file, align_file, cmp_file;
  opts.add ("ps", param_set, "choose parameter set");
  opts.add ("lf", fold_file, "load fold parameters from Telegraph file", false);
  opts.add ("la", align_file, "load align parameters from Telegraph file", false);
  opts.add ("lc", cmp_file, "load cmp parameters from Telegraph file", false);

  try
    {
      Gramset gramset;

      opts.parse_or_die();
      if (opts.args.size())
	{
	  SExpr_file sexpr_file (opts.args[0].c_str());
	  PCFG_builder::init_gramset (sexpr_file.sexpr.find_or_die (CFG_GRAMSET), gramset);
	}
      else
	{
	  Stemloc_gramset stemloc_gramset (false);

	  if (fold_file.size())
	    {
	      Telegraph_PScores_adaptor ss_tgio (stemloc_gramset.single_scfg_pscores);
	      ss_tgio.read (fold_file.c_str(), true);
	    }

	  if (align_file.size())
	    {
	      if (param_set.empty() || stemloc_gramset.phmm_cfg_map.find (param_set) == stemloc_gramset.phmm_cfg_map.end())
		THROWEXPR ("You need to specify a valid param_set");
	      Telegraph_PScores_adaptor qa_tgio (stemloc_gramset.phmm_cfg_map[param_set].hmm_pscores);
	      qa_tgio.read (align_file.c_str(), true);
	    }

	  if (cmp_file.size())
	    {
	      if (param_set.empty() || stemloc_gramset.phmm_cfg_map.find (param_set) == stemloc_gramset.phmm_cfg_map.end())
		THROWEXPR ("You need to specify a valid param_set");
	      Telegraph_PScores_adaptor sp_tgio (stemloc_gramset.phmm_cfg_map[param_set].cfg_pscores);
	      sp_tgio.read (cmp_file.c_str(), true);
	    }

	  stemloc_gramset.fix_group_suffix();

	  gramset = stemloc_gramset;
	}

      Stemloc_gramset::prefix_comments (cout);
      PCFG_builder::gramset2stream (cout, gramset);

    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
