#include "seq/suffix.h"
#include "util/opts_list.h"

int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <sequence file>",
		  "print word counts using a suffix tree");

  int min_word_len;
  int max_word_len;
  int min_count;

  opts.newline();
  opts.add ("l -length", min_word_len = 12, "minimum word length");
  opts.add ("m -maxlen", max_word_len = 100, "maximum word length, 0 to unlimit");
  opts.add ("c -count", min_count = 10, "minimum count");

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      sstring seq_file = opts.args[0];

      // read in sequences
      FASTA_sequence_database seq_db (seq_file.c_str(), 0, Profile_flags_enum::DSQ);
      const Alphabet& alphabet = seq_db.detect_alphabet();

      // make suffix tree
      Suffix_tree tree (seq_db.index, max_word_len, min_word_len, min_count);
      if (CTAGGING(3,SUFFIX)) tree.dump (CL, &alphabet);

      // iterate through nodes
      for_iterator (Suffix_tree::Postorder, iter, tree.begin(), tree.end())
	{
	  if (CTAGGING(2,SUFFIX)) CL << "(visiting node " << iter.node_idx() << ")\n";
	  const int count = (*iter).count;
	  if (iter.suffix_len() >= min_word_len && count >= min_count)
	    {
	      // convert dsq back to text
	      const Digitized_biosequence dsq (iter.suffix());
	      Biosequence seq;
	      alphabet.dsq2seq (dsq, seq);
	      // output
	      cout << seq << " " << count << "\n";
	    }
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
