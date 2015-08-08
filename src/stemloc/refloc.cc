#include "scfg/cfgsexpr.h"
#include "stemloc/quickalign.h"
#include "seq/pkeywords.h"

#define REFLOC_HMM_TAG "hmm"

int main (int argc, char** argv) {
  try {
    INIT_OPTS_LIST (opts, argc, argv, 2,
		    "[options] <refs.fasta> <reads.fasta>",
		    "configurable pairwise alignment");

    opts.newline();
    opts.print_title ("General options");

    sstring model_filename, trained_model_filename;
    bool do_alignment;

    opts.newline();
    opts.add ("m -model", model_filename = "", "read model & parameters from file", false);
    opts.add ("t -train", trained_model_filename = "", "do training, output trained model & parameters to file", false);
    opts.add ("a -align", do_alignment = true, "align sequences & output alignment", true);

    opts.parse_or_die();
    const sstring& refs_filename = opts.args[0];
    const sstring& reads_filename = opts.args[1];

    Alphabet& alphabet = CFG_alphabet;  // hardwired for now
    Pair_PHMM hmm (0, alphabet);
    PScores pscores;
    Alphabet_group null_emit;
    set<int> pad_states;
    set<int> mutable_pgroups;

    Quick_align* quick_align = NULL;

    if (model_filename.size()) {
      SExpr_file model_file (model_filename.c_str());
      hmm = PCFG_builder::init_hmm (model_file.sexpr, pscores, null_emit, pad_states, mutable_pgroups, false);
    } else {
      quick_align = new Quick_align (pscores, null_emit);
      mutable_pgroups = quick_align->mutable_pgroups;
      set<int> non_pad = quick_align->non_pad_states();
      for (int s = 0; s < quick_align->states(); ++s)
	if (non_pad.find(s) == non_pad.end())
	  pad_states.insert (s);
      hmm = *quick_align;
    }

    FASTA_sequence_database refs_db (refs_filename.c_str(), &alphabet);
    FASTA_sequence_database reads_db (reads_filename.c_str(), &alphabet);

    if (trained_model_filename.size()) {
      // do training here

      ofstream out (trained_model_filename.c_str());
      PCFG_builder::hmm2stream (out, hmm, pscores, mutable_pgroups, REFLOC_HMM_TAG, " ", null_emit.group_idx, &pad_states);
    }
  
    if (do_alignment) {
      // do alignment here
      // output alignment
    }

    if (quick_align)
      delete quick_align;
  
  } catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
