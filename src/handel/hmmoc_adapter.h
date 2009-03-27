#ifndef HMMOC_ADAPTER_INCLUDED
#define HMMOC_ADAPTER_INCLUDED

#include "handel/movement.h"
#include "handel/recorder.h"

// default stem for generated C++ header & code filenames
#define HMMOC_CPP_FILENAME_PREFIX "handel2hmmoc"

// filename suffices
#define HMMOC_XML_SUFFIX    ".xml"
#define HMMOC_CC_SUFFIX     ".cc"
#define HMMOC_H_SUFFIX      ".h"
#define HMMOC_PARAMS_SUFFIX ".dat"

// cache of compiled binaries
struct Compiled_executable_cache : map<sstring,sstring>
{
  // data
  bool loaded;  // flag indicating whether cache has already been loaded
  unsigned int loaded_size;  // number of loaded entries
  // destructor
  Compiled_executable_cache() : loaded(false), loaded_size(0) { }
  ~Compiled_executable_cache();
  // load & save methods
  void save_cache() const;
  void load_cache();
};

// Handel-->HMMoC adapter
struct HMMoC_adapter : Grammar_state_enum
{
  // static data
  static sstring root_directory, gcc_exec, gcc_args, tmp_directory, cache_filename;
  static bool leave_debris;
  static int two_d_banding, three_d_banding, four_d_banding;  // diameters of band constraint for 2D, 3D & 4D alignment
  static int model_id;
  static Compiled_executable_cache compiled_executable;

  // member data, summarizing Handel_movement reference
  Handel_movement& move;
  Transducer_SExpr_file& peeled;
  const int alphSize;
  const ETree& etree;
  vector<Score_profile*> peeled_seq;
  long band_diameter;

  // code generation options
  bool generate_forward, generate_backward, generate_sample, generate_viterbi, generate_baum_welch;

  // mapping from acyclic transducer to hmmoc namespace
  vector<int> peeled_obs;
  vector<sstring> peeled_seq_name, peeled_seq_id, peeled_len_id, peeled_sym_id, peeled_index_id, peeled_dummy_id;

  // types for representing transition & emission param indices
  // transitions
  struct Transition_param : pair<int,int>
  {
    // constructors
    Transition_param() { }
    Transition_param (int src, int dest);
    // readout
    Prob eval (HMMoC_adapter& adapter);
  };

  // emissions
  struct Emission_param
  {
    enum Param_type { Match = 0, Insert = 1, Delete = 2 };
    ENode branch_trans_index;
    Param_type param_type;
    int state, sym1, sym2;
    // constructors
    Emission_param() { }
    Emission_param (ENode branch_trans_index,
		    Param_type param_type,
		    int state,
		    int sym1,
		    int sym2 = -1);
    // readout
    Prob eval (HMMoC_adapter& adapter);
  };

  // transition & emission parameter indices
  vector<Transition_param> trans_vec;
  vector<Emission_param> emit_vec;

  // methods

  // constructor
  HMMoC_adapter (Handel_movement& move);

  // dump acylic transducer as a HMMoC model string
  // calling this method sets up trans_vec and emit_vec
  void dump_hmmoc_model (ostream& out,
			 const char* cpp_filename_prefix = HMMOC_CPP_FILENAME_PREFIX,
			 const char* cpp_header_filename_prefix = HMMOC_CPP_FILENAME_PREFIX,  /* if null, all code goes in one .cc file */
			 const char* additional_code = 0);

  // dump acyclic transducer to file & return it as a string
  sstring dump_hmmoc_model_to_file (const char* hmmoc_filename_prefix = HMMOC_CPP_FILENAME_PREFIX,
				    const char* cpp_filename_prefix = HMMOC_CPP_FILENAME_PREFIX,
				    const char* cpp_header_filename_prefix = HMMOC_CPP_FILENAME_PREFIX,
				    const char* additional_code = 0);

  // dump params to file
  void dump_hmmoc_params_to_file (const char* param_filename);

  // helper method to write parameter to parameter stream
  void write_param (Prob param, ostream& param_stream);
  void write_param (long param, ostream& param_stream);

  // run HMMoC, run gcc, run generated program and return output as an array of line strings
  // (uses temporary filenames; should eventually clean up after itself, but spoor is useful for debugging)
  bool exec_program (const char* main_code, vector<sstring>& result);  // returns true if succeeded

  // delegate forward algorithm
  // (In fact, the hmmoc code uses the Backward algorithm, since hmmoc doesn't do Forward traceback; the results should be indistinguishable)
  // NB this subroutine calls exec_recursion_and_traceback(), which messes with the generate_* flags
  bool exec_forward (Score& fwd_sc, vector<vector<int> >& sample_path, int n_samples);  // returns true if succeeded

  // delegate Viterbi algorithm
  // NB this subroutine calls exec_recursion_and_traceback(), which messes with the generate_* flags
  bool exec_Viterbi (Score& viterbi_sc, vector<int>& viterbi_path);  // returns true if succeeded

  // common helper method used by exec_forward and exec_Viterbi.
  // fills a DP matrix, does zero or more tracebacks, extracts results.
  // if sample==true, then Backward algorithms are used; otherwise, Viterbi algorithms are used.
  // NB this subroutine messes with the generate_* flags
  bool exec_recursion_and_traceback (Score& final_sc, vector<vector<int> >& sample_path, int n_samples, bool sample);

  // helper: EHMM accessor
  // currently this defaults to the eliminated EHMM, which may not be quite what you want for Viterbi
  EHMM_transducer_scores& ehmm();

  // helper: mapping between acyclic transducer state indices & hmmoc state identifiers
  static sstring state2id (int state_index);
  static int id2state (const sstring& hmmoc_id);

  // helper: hmmoc state emission identifier
  sstring hmmoc_emission_id (int state);
};


#endif /* HMMOC_ADAPTER_INCLUDED */

