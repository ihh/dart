#include "seq/stockholm.h"
#include "newmat/newmat.h"
#include "util/vector_output.h"
#include <string.h>

#define Stockholm_method_tag 			"SimMeth"
#define Stockholm_average_similarity_tag	"AvSim"
#define Stockholm_autosimilarity_tag		"AutoSim"

bool skip_alignment(ifstream& in);

int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, 2, "[options] (AMA|SPS|PPV|TCS) <stockholm database>", "Alignment database centroid tool");
 
  int nth;
  bool find_auto;
  opts.newline(); 
  opts.add ("n -only-nth", nth = 1, "only use every nth alignment in distance matrix");  
  opts.add ("auto", find_auto = true, "output time-averaged alignment similarity as a function of lag");  
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
     // CLOG(4) <<"Loading alignments\n";
      const sstring filename = opts.args[1];
      
      int runlength = 0;
      ifstream file (filename.c_str());
      if (!file) THROWEXPR (filename << ": File not found");
      CLOG(4) <<"Getting runlength\n";
      while ( file && !file.eof())
      	{
	  if (skip_alignment(file)){ 
	 	 runlength++;
	  }
	}
      file.close();
      CLOG(4) <<"Runlength: "<<runlength<<"\n";
      CLOG(4) <<"Creating "<<(int) runlength/nth<<"^2 distance matrix\n";
      CLOG(4) <<"Space required: "<<sizeof(float) * (runlength/nth) * (runlength/nth) <<" bits\n";
      ArrayLengthSpecifier n((int) runlength/nth);
      SquareMatrix distance_matrix(n);
      distance_matrix = 0.0;
      ifstream file1(filename.c_str());
      for(int i = 1; i <= runlength; i++){
	if ((i % nth) > 0){
	  skip_alignment(file1);
	  continue;
	}
	else {
	  Sequence_database seq_db;
	  Stockholm align1; 
	  align1.read_Stockholm(file1, seq_db);
	  ifstream file2(filename.c_str());
	  for(int j = 1; j <= runlength; j++){
	    if((j % nth) > 0){
	      skip_alignment(file2);
	      continue;
	    }
	    if ( j < i ){
	      skip_alignment(file2);
	      continue;
	    }
	    else {
		Stockholm align2;
		Sequence_database sec_db;
		align2.read_Stockholm(file2, sec_db);
		if(method.compare("AMA") == 0){
		  CLOG(3) << "Preforming comparison "<< i/nth<<", "<<j/nth<<"\n";
		  distance_matrix(i / nth, j / nth) = align1.ama_similarity(align2);
		  distance_matrix(j / nth, i / nth) = distance_matrix(i / nth ,j / nth);
		}
	        else if(method.compare("SPS") == 0){	
		  const int res_pairs12 = align1.residue_pair_overlap (align2);
		  const int res_pairs1 = align1.residue_pairs();
		  const int res_pairs2 = align2.residue_pairs();
		  distance_matrix(i / nth, j / nth) = (float) res_pairs12 / (float) res_pairs1;
		  distance_matrix(j / nth, i / nth) = (float) res_pairs12 / (float) res_pairs2;
		}
	 	else if(method.compare("PPV") == 0){
		  const int res_pairs12 = align1.residue_pair_overlap (align2);
		  const int res_pairs1 = align1.residue_pairs();
		  const int res_pairs2 = align2.residue_pairs();
		  distance_matrix(i / nth, j / nth) = (float) res_pairs12 / (float) res_pairs2;
		  distance_matrix(j / nth, i / nth) = (float) res_pairs12 / (float) res_pairs1;
	  	}
	  	else if(method.compare("TCS") == 0){
	   	  const int cols1 = align1.columns();
		  const int cols2 = align2.columns();
		  const int col_overlap = align1.column_overlap(align2);
		  distance_matrix(i / nth, j / nth) = col_overlap/cols1;
		  distance_matrix(j / nth, i / nth) = col_overlap/cols2;
		}
	  	else { //sanity check, should never reach here
	   	  THROWEXPR ("Invalid comparison method!");
	  	}
             }
	  }
	  file2.close();
	}
     } 
     file1.close();    
     // find centroid
     float max = 0.0;
     int max_index = -1;
     CLOG(4) << "Finding centroid...\n";
     for(int x = 1; x <= runlength / nth; x++){
     	const float this_sum = distance_matrix.Row(x).Sum();
	if (max_index < 0 || this_sum > max){
		max = this_sum;
		max_index = x;
	}
     }
     max = max / runlength;

     // create alignment; annotate similarity & avg similarity
     Stockholm final;
     ifstream final_file(filename.c_str());
     int index = 1;
     while ( final_file && !final_file.eof() && (index < max_index)){
	skip_alignment(final_file);
	index++;
     }
     Sequence_database new_db;
     final.read_Stockholm(final_file, new_db);
     sstring tmp;
     final.add_gf_annot(Stockholm_method_tag, method);  
     tmp << max;
     final.add_gf_annot(Stockholm_average_similarity_tag, tmp);  

     // annotate alignment with stationary autocorrelation-like measure, if requested
     if (find_auto) {
       vector<double> autovec;
       const int matrix_size = runlength / nth;
       for (int lag = 0; lag < matrix_size; lag++) {
	 double total = 0.;
	 for (int i = 1; i <= matrix_size - lag; i++)
	   total += distance_matrix (i, i + lag);
	 const double time_avg = total / (double) (matrix_size - lag);
	 autovec.push_back (time_avg);
       }
       sstring autovec_string;
       autovec_string << autovec;
       final.add_gf_annot(Stockholm_autosimilarity_tag, autovec_string);
     }

      // output final alignment
     final.write_Stockholm(cout); 
 
  }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }
  return 0;
}

bool skip_alignment (ifstream& in){
  Regexp re_sep = "^[ \t]*//[ \t]*$"; 
  sstring input;
  while(!re_sep.Match(input.c_str())){
	if (in.eof()){
		return false;
	}
	input.getline(in).chomp();
  }
  return true;
}
