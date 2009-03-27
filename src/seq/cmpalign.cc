#include "seq/stockholm.h"

int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, 2, "[options] <reference alignment> <comparison alignment>", "Alignment comparison tool");
  
  //bool flat_format;
  bool scores_only;

  opts.newline();
  //opts.add ("f -flat-format", flat_format = false, "output alignment statistics in a flat-format");
  opts.add ("s -scores-only", scores_only = false, "output scores only, in AMA\\tSPS\\tPPV\\tTCS format");

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
    //  opts.parse_or_die();
      
      CLOG(4) <<"Loading alignments\n";
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
      align1.erase_empty_columns();
      align2.erase_empty_columns();
      const int res_pairs1 = align1.residue_pairs();
      const int res_pairs2 = align2.residue_pairs();
      const int in_pairs1 = align1.count_pairwise_insertions();
      const int del_pairs1 = align1.count_pairwise_deletions();
      const int in_pairs2 = align2.count_pairwise_insertions();
      const int del_pairs2 = align2.count_pairwise_deletions();
      const int in_pairs12 = align1.pairwise_insertion_overlap (align2);
      const int del_pairs12 = align1.pairwise_deletion_overlap (align2);
      const int res_pairs12 = align1.residue_pair_overlap (align2);
      const int cols1 = align1.columns();
      const int cols2 = align2.columns();
      const int col_overlap = align1.column_overlap(align2);
      const float sps_score = (float) res_pairs12 / (float) res_pairs1;
      const float ppv_score = (float) res_pairs12 / (float) res_pairs2;
      const float ama_score = align1.ama_similarity(align2);
      const float tc_score = (float) col_overlap / (float) cols1;
	
	if (scores_only){
		cout << ama_score <<"\t"<<sps_score<<"\t"<<ppv_score<<"\t"<<tc_score<<"\n";
		return 0;
	}
			

      cout << "Alignment '" << argv[1] << "' has "<< cols1 <<" columns, " << res_pairs1 << " residue pairs" << " and " << in_pairs1 << " insertions and "<<del_pairs1<<" deletions\n";
      cout << "Alignment '" << argv[2] << "' has "<< cols2 <<" columns, "<< res_pairs2 << " residue pairs" << " and " << in_pairs2 << " insertions and "<<del_pairs2<<" deletions\n";
      cout << "Overlap: " << col_overlap << " columns, "<<  res_pairs12 << " residue pairs " << in_pairs12 << " insertions and "<<del_pairs12<< " deletions\n";
      cout << "Normalized AMA Similarity Score:\t"<< ama_score << "\n";
      cout << "Sum of Pairs Score:             \t"<< sps_score << "\n";
      cout << "Positive Predictive Value:      \t"<< ppv_score << "\n";
      cout << "Total Column Score:             \t"<< tc_score << "\n";
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }
  return 0;
}

