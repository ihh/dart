#include "seq/stockholm.h"
#include "newmat/newmat.h"
#include <string.h>

int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, 2, "[options] (AMA|SPS|PPV|TCS) <stockholm database>", "Alignment database centroid tool");
 
  int nth = 1;
  opts.newline(); 
  opts.add ("n -only-nth", nth = 1, "only use every nth alignment in distance matrix");  
  //bool flat_format;
  //bool scores_only;

  //opts.newline();
  //opts.add ("f -flat-format", flat_format = false, "output alignment statistics in a flat-format");
  //opts.add ("s -scores-only", scores_only = false, "output scores only, in AMA\\tSPS\\tPPV\\tTCS format");

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
      const string method = opts.args[0];
      if(!((method.compare("AMA")&&method.compare("SPS")&&method.compare("PPV")&&method.compare("TCS")) == 0)){
	CLOGERR << "ERROR: Invalid Comparison Measure: "<<method<<"\n"<< opts.short_help(); exit(1);
      } 
      CLOG(4) <<"Using comparison measure: "<<method<<"\n";
	
      CLOG(4) <<"Loading alignments\n";
      const sstring filename = opts.args[1];
      
      Stockholm_database stock_db, tmp_db; 
      Stockholm align1, align2;
      Sequence_database seq_db;

      ifstream file (filename.c_str());
      if (!file) THROWEXPR (filename << ": File not found");
      stock_db.read (file, seq_db);
      file.close();
      tmp_db = stock_db;
     
      CLOG(4) <<"Creating "<<(int) stock_db.size()/nth<<"^2 distance matrix\n";
      ArrayLengthSpecifier n(stock_db.size());
      SquareMatrix distance_matrix(n);
      int i = 1;
      int j = 1;
      for_const_contents(list<Stockholm>, stock_db.align, align1_tmp){
	if ((i % nth) > 0) {
		i++;
		if(i > stock_db.align.size()){
		  break;
		}
		continue;
	}
      	CLOG(3) <<"Computing row "<<i / nth<<"\n";
	distance_matrix(i / nth, i / nth ) = 1;
      	tmp_db.align.pop_front();
	j = i + 1;
      	for_const_contents(list<Stockholm>, tmp_db.align, align2_tmp){
	   if((j % nth) > 0) {
		j++;
		if(j > stock_db.align.size()){
		  break;
		}
		continue;
	   }
	   Stockholm align1 (*align1_tmp);
	   Stockholm align2 (*align2_tmp);
      	   align1.erase_empty_columns();
      	   align2.erase_empty_columns();
	   CLOG(1)<< "Computing matrix entries: "<<i / nth<<","<<j / nth<<"\n";
	   if(method.compare("AMA") == 0){
		distance_matrix(i / nth, j / nth) = align1.ama_similarity(align2);
		distance_matrix(j / nth, i / nth) = distance_matrix(i,j);
	   }
	   else if(method.compare("SPS") == 0){
	   	const int res_pairs12 = align1.residue_pair_overlap (align2);
		const int res_pairs1 = align1.residue_pairs();
		const int res_pairs2 = align2.residue_pairs();
    	   	distance_matrix(i, j) = (float) res_pairs12 / (float) res_pairs1;
		distance_matrix(j, i) = (float) res_pairs12 / (float) res_pairs2;
	   }
	   else if(method.compare("PPV") == 0){
	   	const int res_pairs12 = align1.residue_pair_overlap (align2);
		const int res_pairs1 = align1.residue_pairs();
		const int res_pairs2 = align2.residue_pairs();
    	   	distance_matrix(i, j) = (float) res_pairs12 / (float) res_pairs2;
		distance_matrix(j, i) = (float) res_pairs12 / (float) res_pairs1;
	   }
	   else if(method.compare("TCS") == 0){
	   	const int cols1 = align1.columns();
           	const int cols2 = align2.columns();
           	const int col_overlap = align1.column_overlap(align2);
		distance_matrix(i,j) = col_overlap/cols1;
		distance_matrix(j,i) = col_overlap/cols2;
	   }
	   else { //sanity check, should never reach here
	   	THROWEXPR ("Invalid comparison method!");
	   }
	   j++;
	}
	i++;
     }
     
     float max = 0.0;
     int max_index;
     CLOG(4) << "Finding centroid...\n";
     for(int x = 1; x <= stock_db.size() / nth; x++){
     	const float this_sum = distance_matrix.Row(x).Sum();
	if (this_sum > max){
		max = this_sum;
		max_index = x;
	}
     }
     max = max / stock_db.size();
     Stockholm final (*stock_db.align_index[max_index-1]);
     sstring tmp;
     final.add_gf_annot("Similarity Measure:", method);
     tmp << max;
     final.add_gf_annot("Average Similarity:", tmp);
     final.write_Stockholm(cout); 
 
  }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }
  return 0;
}

