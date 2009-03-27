// pk 3/05 from xrate.cc support asymmetric matrices
// pk 3/05 accept long input lines in stockholm file
// pk 4/05 eliminate the "--quick" and "--rind" command line options
// pk 4/05 consolidate

#include "irrev/irrev_em_matrix.h"
#include "tree/subdistmat.h"
#include "tree/pam.h"
#include "tree/hasegawa.h"

// transition/transversion ratio for HKY model
#define HKY_IV_RATIO 10.

// no. of processes to fork (set to 1 to disable this feature, as pipes seem to be broken under gcc 3.3)
#define max_fork 1

// main program
int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <alignment database in Stockholm format>",
		  "estimate a substitution rate matrix from an alignment database using Expectation Maximisation\n");

  int C;
  bool train;
  bool intra;
  bool inter;
  bool rndinit;
  bool infer_classes;
  double tres;
  //  int max_fork;
  int max_rounds;
  double em_min_inc;
  int forgive;
  bool norm;
  sstring alph_chars;
  sstring init_hsm_filename;

  opts.newline();
  opts.newline();
  opts.print ("General xrate options\n");
  opts.print ("---------------------\n");
  opts.add ("i init -initialise", init_hsm_filename = "", "read initial rate matrix from this file", FALSE);
  opts.add ("ri -rndinit", rndinit = FALSE, "randomise matrix if no initialisation file supplied");
  opts.add ("a alph -alphabet", alph_chars = "", "use alphabet with these characters", FALSE);
  opts.add ("t train -train", train = 1, "train matrices on indexed alignments using EM");
  opts.add ("tr tres -timeres", tres = .01, "fractional resolution of branch lengths");
  opts.add ("n -normalise", norm = FALSE, "normalise rate matrix");
  opts.add ("hsmhelp", &EM_matrix_base::hsm_help, "\thelp on rate matrix file format");

  opts.newline();
  opts.newline();
  opts.print ("EM options\n");
  opts.print ("----------\n");
  opts.add ("g forgive -forgive", forgive = 0, "number of bad rounds of EM to forgive");
  opts.add ("mi -mininc", em_min_inc = .001, "minimum fractional increase in log-likelihood per round of EM");
  opts.add ("mr -maxrounds", max_rounds = -1, "max number of rounds of EM, or -1 to unlimit");
  //  opts.add ("f fork -fork", max_fork = 1, "number of processes to fork [KNOWN BUG: avoid using this feature, as pipes appear to hang under gcc 3.3.1]");

  opts.newline();
  opts.newline();
  opts.print ("Hidden site class options\n");
  opts.print ("-------------------------\n");
  opts.add ("c classes -classes", C = 1, "number of dynamic hidden site classes");
  opts.add ("intra -intra", intra = 1, "train inter-class substitution rates");
  opts.add ("inter -inter", inter = 1, "train inter-class substitution rates");
  opts.add ("ic -inferclass", infer_classes = FALSE, "infer class labels");

  opts.newline();
  opts.newline();
  opts.print ("Acknowledgements\n");
  opts.print ("----------------\n");
  opts.print ("Released under the GNU Public License by Ian Holmes, Department of Bioengineering, University of California, Berkeley. Email <ihh@berkeley.edu>\n");
  opts.print ("Programmers: Ian Holmes, Pete Klosterman, Gerton Lunter, Robert Davies.\n");
  opts.print ("Beta testers: Carolin Kosiol, Dawn Brooks, Caleb Webber, Nick Goldman, Simon Whelan.\n");
  opts.newline();
  opts.print ("For information about Stockholm alignment format, see http://www.cgr.ki.se/cgb/groups/sonnhammer/Stockholm.html\n");
  opts.newline();
  opts.print ("Reference\n");
  opts.print ("---------\n");
  opts.print ("I.Holmes and G.M.Rubin, 2002. An Expectation Maximization algorithm for training hidden substitution models. Journal of Molecular Biology 317:5, 757-768.\n");
  opts.newline();

  // parse the command line
  //
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  try
    {
      // get args
      const char* alignment_db_filename = opts.args[0].c_str();

      // initialise database
      Sequence_database seq_db;
      Stockholm_database stock_db;
      ifstream align_db_in (alignment_db_filename);
      if (!align_db_in) THROWEXPR ("Couldn't open alignment file '" << alignment_db_filename << "'");
      stock_db.read (align_db_in, seq_db, 12345678);			// pk prev. max line length was 123456 set in sstring.h
      Tree_alignment_database align_db (seq_db);
      align_db.initialise_from_Stockholm_database (stock_db, FALSE);

      // create dummy alphabet
      Alphabet dummy_alph ("dummy_alphabet", alph_chars.size());
      dummy_alph.init_chars (alph_chars.c_str());

      // get the initial alphabet
      const Alphabet& init_alph = dummy_alph.size() ? dummy_alph : seq_db.detect_alphabet();

      // initialise EM_matrix size, alphabet & number of rounds
      // (the effort of initialising the size & alphabet will be wasted if we later read an initial matrix from a file, but never mind)
      Irrev_EM_matrix hsm (C, init_alph.size(), max_fork);
      hsm.init_alphabet (init_alph);
      hsm.em_min_inc = em_min_inc;
      if (max_rounds >= 0)
	hsm.em_max_iter = max_rounds;

      // read in initial matrix, if specified; otherwise, randomise
      bool initially_randomised = FALSE;
      if (init_hsm_filename.size())
	  {
		ifstream hsm_file (init_hsm_filename.c_str());
		if (!hsm_file) THROWEXPR ("HSM file '" << init_hsm_filename << "' not found");
		hsm.read (hsm_file);
	  }
      else
	  {
		// some heuristic rules for the randomness of the initial randomised matrix
		// (chosen so that expected total substitution rate per site is 1)
		const double prior_dev = rndinit ? (.5 / (double) (C*init_alph.size())) : 0.;
		const double inter_intra_ratio = .1;
		const double intra_min = 1. / (double) (inter_intra_ratio * (double) (C-1) + (init_alph.size()-1));
		const double intra_dev = rndinit ? (intra_min / 5.) : 0.;
		const double inter_min = inter_intra_ratio * intra_min;
		const double inter_dev = rndinit ? (inter_intra_ratio * intra_dev) : 0.;
		hsm.randomise (prior_dev, intra_min, intra_dev, inter_min, inter_dev);
		initially_randomised = TRUE;
	  }

      // estimate missing trees in alignment database
      // firstly, get a substitution model
      Substitution_matrix_factory* submat_factory = 0;
      if (initially_randomised)  // if EM_matrix is randomised, then create a temporary "realistic" substitution model
	  {
		// if this is a known alphabet, use a plausible matrix (PAM for protein, HKY for DNA/RNA)
		if (&init_alph == &Protein_alphabet)
		  submat_factory = new PAM_factory;
		else if (&init_alph == &DNA_alphabet || &init_alph == &RNA_alphabet)
	    {
	      seq_db.seqs2dsqs (init_alph);  // digitize the sequences
	      const double base_pseudocount = 1.;
	      const vector<Prob> base_composition = seq_db.get_null_model (init_alph.size(), base_pseudocount);
	      seq_db.clear_dsqs();  // probably not necessary, but ensures wrong alphabet is not used later on
	      submat_factory = new Hasegawa_matrix_factory (1. / (1. + HKY_IV_RATIO),
							    HKY_IV_RATIO / (1. + HKY_IV_RATIO),
							    base_composition, TRUE);
	    }
		else // otherwise, use a "flat" matrix
	    {
	      Irrev_EM_matrix* flat_em_mat = new Irrev_EM_matrix (1, init_alph.size(), max_fork);
	      flat_em_mat->randomise (0., 1./(double)init_alph.size(), 0., 0., 0.);
	      flat_em_mat->init_alphabet (init_alph);
	      submat_factory = flat_em_mat;
	    }
	  }
      else
		submat_factory = &hsm;

      // secondly, estimate the trees
      seq_db.seqs2scores (submat_factory->alphabet());  // convert sequences to Score_profile's
      Subst_dist_func_factory dist_func_factory (*submat_factory);
      align_db.estimate_missing_trees_by_nj (dist_func_factory);
      seq_db.clear_scores();  // clear the Score_profile's to ensure wrong alphabet not used later

      // finally, delete the temporary Substitution_matrix_factory if appropriate
      if (initially_randomised)
		delete submat_factory;

      // now we can tell the EM_matrix about the Tree_alignment_database
      hsm.init_cache (&align_db, tres);
      hsm.update();

      // digitise sequences
      const Alphabet& hidden_alphabet = hsm.alphabet();
      seq_db.seqs2scores (hidden_alphabet);

      // train matrix
      Loge log_likelihood;
      if (train)
		log_likelihood = hsm.iterate_quick_EM (intra, inter, forgive, FALSE);	// rind not an option
      else
		log_likelihood = hsm.log_likelihood();

      // display class labels
      if (infer_classes)
		hsm.display_classes (cout);

      // show log-likelihood
      CTAG(7,XRATE) << "Final log-likelihood = " << Nats2Bits(log_likelihood) << " bits\n";

      // normalise rate matrix
      if (norm)
	  {
		CTAG(7,XRATE) << "Normalising rates to 1 expected substitution per site per unit time\n";
		hsm.normalise_substitution_rate();
	  }

      // output matrix
      hsm.write (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << "ERROR: " << e.what();
      exit(1);
    }

  return 0;
}
