#include "seq/stockholm.h"
#include "util/vector_output.h"

typedef pair<int,int> Base_pair;
typedef set<Base_pair> Base_pair_set;

Base_pair_set fold_string_to_base_pair_set (const sstring& fold_string)
{
  Base_pair_set pairs;
  stack<int> lpos;
  for (int pos = 0; pos < (int) fold_string.size(); ++pos)
    {
      const char c = fold_string[pos];
      if (Fold_char_enum::is_lchar(c))
	lpos.push (c);
      else if (Fold_char_enum::is_rchar(c))
	{
	  if (lpos.empty())
	    {
	      THROWEXPR ("Too many >'s");
	    }
	  else
	    {
	      pairs.insert (Base_pair (lpos.top(), pos));
	      lpos.pop();
	    }
	}
    }
  if (lpos.size())
    THROWEXPR ("Too many <'s");
  return pairs;
}

int main (int argc, char** argv)
{
  INIT_OPTS_LIST(opts,argc,argv,2,"[options] <alignment#1> <alignment#2>","RNA structural alignment comparison tool");

  try
    {
      // parse command line
      opts.parse_or_die();

      // load alignments
      CLOG(4) << "Loading alignments\n";
      const sstring filename1 = opts.args[0];
      const sstring filename2 = opts.args[1];

      Stockholm align1, align2;
      Sequence_database seq_db;

      ifstream file1 (filename1.c_str());
      if (!file1) THROWEXPR (filename1 << ": File not found");
      align1.read_Stockholm (file1, seq_db);
      file1.close();

      ifstream file2 (filename2.c_str());
      if (!file2) THROWEXPR (filename2 << ": File not found");
      align2.read_Stockholm (file2, seq_db);
      file2.close();

      // put all results in a single vector
      vector<int> results;

      // get aligned-pair overlap
      CLOG(5) << "Checking alignment overlap\n";
      const int apairs1 = align1.residue_pairs();
      const int apairs2 = align2.residue_pairs();
      const int apairs12 = align1.residue_pair_overlap (align2);

      results.push_back (apairs1);
      results.push_back (apairs2);
      results.push_back (apairs12);

      // output results
      cout << "Residue pairs: ";
      cout << apairs1 << " ('" << filename1 << "'), ";
      cout << apairs2 << " ('" << filename2 << "'), ";
      cout << "overlap " << apairs12 << "\n";

      // get sequence names
      set<sstring> seqname;
      for_const_contents (vector<sstring>, align1.row_name, name)
	seqname.insert (*name);
      for_const_contents (vector<sstring>, align2.row_name, name)
	seqname.insert (*name);

      // get base-pair sets and overlap, and output those
      CLOG(5) << "Checking structure overlap\n";
      for_const_contents (set<sstring>, seqname, name)
	{
	  const sstring fold1 = align1.get_fold_string (*name);
	  const sstring fold2 = align2.get_fold_string (*name);

	  Base_pair_set pairs1 = fold_string_to_base_pair_set (fold1);
	  Base_pair_set pairs2 = fold_string_to_base_pair_set (fold2);

	  int overlap = 0;
	  for_const_contents (Base_pair_set, pairs1, bp)
	    if (pairs2.find(*bp) != pairs2.end())
	      ++overlap;

	  cout << "Basepairs (sequence '" << *name << "'): ";
	  cout << pairs1.size() << " ('" << filename1 << "'), ";
	  cout << pairs2.size() << " ('" << filename2 << "'), ";
	  cout << "overlap " << overlap << "\n";

	  results.push_back (pairs1.size());
	  results.push_back (pairs2.size());
	  results.push_back (overlap);
	}

      // summary
      cout << "Summary: " << results << '\n';
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }
  return 0;
}

