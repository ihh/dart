#include <iomanip.h>
#include "kimono/trues.h"
#include "util/logfile.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  Opts_list opts (argc, argv, CLOGERR);
  opts.short_description = "Gibbs sampler";
  opts.syntax = "[options] <sequence file>";
  Rnd::add_opts (opts);
  opts.newline();
  Log_stream::add_opts (opts);

  int     motif_len;
  double  motif_del_prob;
  bool    allow_revcomp;
  double  pseudocount_mul;
  int     max_rounds;
  sstring best_file;
  sstring cluster_file;
  sstring initgff;
  bool    stutter;

  opts.newline();
  opts.add ("motiflen",   motif_len = 10,             "motif length");
  opts.add ("revcomp",    allow_revcomp = 1,          "\tallow reverse strand alignments for DNA");
  opts.add ("mdelprob",   motif_del_prob = .0001,     "probability of deleting a motif residue");
  opts.add ("stutter",    stutter = 1,                "use stuttered null model");
  opts.add ("pseudomul",  pseudocount_mul = .1,       "multiplier to convert null model into alignment pseudocounts");
  opts.add ("rounds",     max_rounds = 100,           "\tnumber of rounds of sampling algorithm");
  opts.add ("cfile",      cluster_file = "",          "\tclustering file (one cluster per line: 'CLUSTERNAME MEMBER1 MEMBER2 ...')", 0);
  opts.add ("bestfile",   best_file = "",             "file to write best intermediate results into", 0);
  opts.add ("initgff",    initgff = "",               "GFF file containing initial motif co-ords", 0);
  
  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 1) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // get parameters

      sstring seq_file = opts.args[0];

      // read in the sequences; detect alphabet; digitise; revcomp if allowed

      FASTA_sequence_database seq_db (seq_file.c_str());
      const Alphabet& alphabet = seq_db.detect_alphabet();
      seq_db.seqs2dsqs (alphabet);

      Sequence_database revcomp_db;
      if (alphabet.has_complement() && allow_revcomp) revcomp_db = seq_db.revcomp (alphabet);

      // set up database indices

      Sequence_database_index db_index (seq_db);
      Sequence_database_index revcomp_index (revcomp_db);

      // read in the clusters
      
      vector<vector<sstring> > members;
      vector<sstring> cluster_name;
      vector<int> which_cluster (seq_db.size(), (int) 0);
      if (cluster_file != "") {
	ifstream cf (cluster_file.c_str());
	sstring c;
	while (1) {
	  c.getline (cf);
	  vector<sstring> f = c.split();
	  if (f.size() == 0) break;
	  cluster_name.push_back (f[0]);
	  f.erase (f.begin());
	  members.push_back (f);
	}
      } else {
	members.push_back (vector<sstring> ());
	for_const_contents (Sequence_database, seq_db, s) members.back().push_back ((*s).name);
	cluster_name.push_back ("<default-cluster>");
      }
      
      int k = cluster_name.size();

      // do each cluster successively
      
      for (int i = 0; i < k; ++i) {
	const vector<sstring>& cluster = members[i];

	Gibbs_dataset     gibbs_dataset (db_index, revcomp_index, cluster, alphabet);
	vector<double>    compositional_pseudocounts (alphabet.size(), pseudocount_mul);    // these counts are added on to the composition vector in the null model
	Gibbs_null_model  gibbs_null (gibbs_dataset, alphabet, compositional_pseudocounts, motif_del_prob, stutter);
	vector<double>    model_pseudocounts = gibbs_null.composition;                      // these are the actual pseudocounts, derived from the composition vector
	for_contents (vector<double>, model_pseudocounts, pc) *pc = *pc * pseudocount_mul;
	Gibbs_alignment_factory  gibbs_factory (gibbs_dataset, alphabet, gibbs_null, model_pseudocounts, motif_len);

	Gibbs_alignment* gibbs = gibbs_factory.new_model();
	for (int g = 0; g < cluster.size(); ++g) gibbs->set_membership_probability (g, 0);

	// best alignment to date:
	
	Gibbs_alignment* best = 0;
	int best_sc = -InfinityScore;

	// read the initialisation GFF file, if there is one

	if (initgff != "") {
	  ifstream initgff_s (initgff.c_str());
	  sstring g;
	  int n = 0;
	  while (!initgff_s.eof()) {
	    g.getline (initgff_s);
	    g.chomp();
	    ++n;
	    vector<sstring> f = g.split("\t");
	    if (f.size() == 0) continue;
	    if (f.size() < 9) { CLOGERR << "Line " << n << " of '" << initgff << "' is not strict GFF; skipping\n"; continue; }
	    int row = gibbs_dataset.row_index (f[0]);
	    if (row >= 0) {
	      bool rev = f[6] == "-";
	      int start;
	      if (rev)
		start = gibbs_dataset.profile[row]->size() - atoi (f[4].c_str());
	      else
		start = atoi (f[3].c_str()) - 1;
	      if (CLOGGING(5))
		CL << "initgff: sequence '" << f[0] << "', " << (rev ? "reverse" : "forward") << " strand, offset " << start << ", motif '" << f[2] << "'\n";
	      gibbs->align_row (row, -start, rev);
	    }
	  }
	}

	CLOG(7) << "Ready to start.\n";

	// start the sampling loop
      
	for (int round = 0; round < max_rounds; ++round)
	  {
	    CLOG(6) << "Round " << round+1 << ": sampling sequence alignments\n";
	    
	    vector<int> sequence_order (cluster.size());
	    for (int g = 0; g < cluster.size(); ++g) sequence_order[g] = g;
	    for (int j = 0; j < cluster.size(); ++j)
	      {
		swap (sequence_order[j], sequence_order [Rnd::rnd_int (sequence_order.size() - j) + j]);
		const int g = sequence_order[j];
		gibbs->sample_row (g);
		if (CLOGGING(3)) gibbs->display (CL);
		if (round > 0) {  // make sure all rows aligned
		  int sc = 0;
		  for (int k = 0; k < cluster.size(); ++k) sc += gibbs->log_likelihood_sc (k);
		  if (sc > best_sc) {
		    if (best) delete best;
		    best = new Gibbs_alignment (*gibbs);
		    best_sc = sc;
		  }
		}
	      }
	    gibbs->sample_motif();
	    gibbs->flush_left();
	  }

	cout << "Cluster '" << cluster_name[i] << "'\n";
	if (best) {
	  best->display (cout);
	  delete best;
	}
	delete gibbs;

	// end of sampling loop
      }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}

