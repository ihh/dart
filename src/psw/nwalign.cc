#include "hmm/pairhmm.h"
#include "seq/stockholm.h"

struct Protein_Needleman_Wunsch_HMM : Pair_HMM_scores
{
  // states, default scores
  enum { Match = 0, Delete = 1, Insert = 2, TotalStates = 3 };
  enum { DefaultMatchScore = +5, DefaultMismatchScore = -4 };

  // method to read BLAST matrix into Match score array
  void read_submat (const char* filename)
  {
    ifstream in (filename);
    if (!in)
      THROWEXPR ("Error opening BLAST substitution matrix file '" << filename << "'");
    read_BLAST_matrix (in, Match);
  }

  // constructor
  Protein_Needleman_Wunsch_HMM (Score gap_open, Score gap_extend)
    : Pair_HMM_scores (TotalStates, &Protein_alphabet)
  {
    // initialise state types
    init_emit (Match, EmitXY, DefaultMismatchScore);
    init_emit (Delete, EmitX, 0);
    init_emit (Insert, EmitY, 0);
    
    // set up the default substitution matrix
    for (int i = 0; i < (int) alphabet->size(); ++i)
      pair_emit[Match](i,i) = DefaultMatchScore;

    // initialise transition scores
    transition (Start, Match) = 0;
    transition (Start, Delete) = -gap_open;
    transition (Start, Insert) = -gap_open;
    transition (Start, End) = 0;

    transition (Match, Match) = 0;
    transition (Match, Delete) = -gap_open;
    transition (Match, Insert) = -gap_open;
    transition (Match, End) = 0;

    transition (Delete, Match) = 0;
    transition (Delete, Delete) = -gap_extend;
    transition (Delete, Insert) = -gap_open;
    transition (Delete, End) = 0;

    transition (Insert, Match) = 0;
    transition (Insert, Insert) = -gap_extend;
    transition (Insert, End) = 0;
  }
};

// main function
int main (int argc, char** argv)
{
  // parse command-line options
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <sequence file>",
		  "protein Needleman-Wunsch aligner");

  sstring gap_penalty;
  Score gap_open;
  Score gap_extend;
  sstring submat_filename;

  opts.print ("Needleman-Wunsch algorithm options\n");
  opts.print ("----------------------------------\n");

  opts.add ("g", gap_penalty = "", "use linear gap penalty", false);
  opts.add ("o", gap_open = 10, "affine gap opening penalty");
  opts.add ("x", gap_extend = 10, "affine gap extension penalty");
  opts.add ("s", submat_filename = "", "read BLAST-format substitution matrix from file", 0);

  opts.parse_or_die();

  // do the work
  try
    {
      // get sequence filename
      const sstring seq_filename = opts.args[0];

      // use linear gap penalty, if specified
      if (gap_penalty.size())
	gap_open = gap_extend = gap_penalty.to_int();

      // make Pair HMM
      Protein_Needleman_Wunsch_HMM nwhmm (gap_open, gap_extend);

      // read in submat
      if (submat_filename.size())
	nwhmm.read_submat (submat_filename.c_str());

      // read in sequences
      FASTA_sequence_database seq_db (seq_filename.c_str(), nwhmm.alphabet);

      // do alignments
      for (int i = 0; i < seq_db.size() - 1; ++i)
	for (int j = i + 1; j < seq_db.size(); ++j)
	  {
	    // get sequences
	    const Named_profile& iseq (seq_db.get_seq (i));
	    const Named_profile& jseq (seq_db.get_seq (j));
	    // make DP matrix
	    const Pair_Viterbi_DP_matrix viterbi_matrix (nwhmm, iseq.prof_sc, jseq.prof_sc);
	    // get alignment
	    const vector<int> viterbi_state_path = viterbi_matrix.optimal_state_path();
	    const Pairwise_path viterbi_alignment_path = nwhmm.convert_state_path_to_alignment (viterbi_state_path);
	    const Alignment viterbi_alignment (viterbi_alignment_path, iseq, jseq);
	    // output
	    viterbi_alignment.write_MUL (cout, *nwhmm.alphabet);
	    cout << Stockholm_file_annotation << ' ' << Stockholm_bit_score_tag << ' ' << viterbi_matrix.final_score << '\n';
	    cout << Stockholm_alignment_separator << '\n';
	  }
    }
  catch (const Dart_exception& e)  // exception; bail out gracefully
    {
      CLOGERR << e.what();
      exit(1);
    }

}
