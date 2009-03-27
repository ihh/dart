#ifndef MATRIX_ARRAY
#define MATRIX_ARRAY

#include "hsm/em_matrix.h"

// Matrix_array is an array of independent EM_matrix's
// with a fixed hidden clique class, conditioned on a hidden sequence class
struct Matrix_array_params
{
  // dimensions
  const Tree_alignment_database* align_db;  // alignment database
  double timepoint_res;  // resolution of time for Timepoint_cache
  int S;  // no. of hidden sequence classes
  int F;  // no. of *fixed* hidden clique classes
  int C;  // no. of *variable* hidden clique classes
  int A;  // alphabet size

  // parameters
  vector<EM_matrix> hsm;
  vector<Loge> sprior;  // log-prior for sequence class
  vector<vector<Loge> > fprior;  // log-prior for fixed clique class, conditioned on sequence class
};

struct Matrix_array : Matrix_array_params
{
  int max_fork;  // this is only stored so EM_matrix can be initialised properly
  // accessors
  inline int m() const { return C*A; }
  const Alphabet& alphabet() { return hsm.front().alphabet(); }

  // constructors/initialisers
  Matrix_array (int S, int F, int C, int A, int max_fork, const Tree_alignment_database* align_db = 0, double timepoint_res = DEFAULT_TIMEPOINT_RES);
  void init_array (int S, int F, int C, int A, int max_fork, const Tree_alignment_database* align_db = 0, double timepoint_res = DEFAULT_TIMEPOINT_RES);
  void init_alphabet (const Alphabet& base_alphabet);
  void update();  // calls update() for all HSM's

  // update statistics
  struct Update_statistics
  {
    const int S;  // no. of hidden sequence classes
    const int F;  // no. of fixed hidden clique classes
    const int states;  // C*A
    Loge log_likelihood;
    vector<EM_matrix::Update_statistics> hsm_stats;  // expected counts for each EM_matrix
    vector<double> scount;  // expected counts for sequence class
    vector<vector<double> > fcount;  // expected counts for fixed clique class
    // methods
    Update_statistics (int S, int F, int states);  // constructor
    void clear();  // zeroes all counts
    void clear_DJU();  // clears DJU's only
    void transform (const Matrix_array& arr, bool symmetrise);  // calls transform() for each hsm_stats
    friend ostream& operator<< (ostream& o, const Update_statistics& stats);
  };
  void equilibrate (Update_statistics& stats) const;
  
  // up/down matrices
  struct Alignment_matrix
  {
    // data
    const Matrix_array& arr;
    const Alignment& align;
    const int S;
    const int F;
    vector<EM_matrix::Alignment_matrix> hsm_matrix;  // Alignment_matrix for each EM_matrix
    Loge total_log_likelihood;  // summed over sequence classes, columns, cliques & clique classes
    // counts, held here as a courtesy to Update_statistics
    vector<double> scount;
    vector<vector<double> > fcount;
    // constructor
    Alignment_matrix (const Matrix_array& arr, int n_align);
    // accessors
    int n_cliques (int column) const;
    // up/down algorithm
    void fill_up();
    void calc_post();  // calculates total_log_likelihood as well as fixed-clique-class posteriors
    void fill_down (Update_statistics& stats);
  };

  // up/down algorithm
  void up_down (Update_statistics& stats, bool likelihood_only, bool symmetrise) const;  // does up-down algorithm on each column
  Update_statistics get_stats_unequilibrated (bool symmetrise) const;
  Update_statistics get_stats() const;
  Loge log_likelihood() const;

  // init method
  void randomise (double seq_dev = .1, double fixed_dev = .1, double prior_dev = .1, double intra_min = .05, double intra_dev = .01, double inter_min = .005, double inter_dev = .001);

  // quick EM
  Update_statistics single_quick_EM (bool rind = 0, bool intra = 1, bool inter = 1);  // does partial M-step only
  Loge iterate_quick_EM (bool rind = 0, bool intra = 1, bool inter = 1, int forgive = 20);  // returns best log-likelihood

  // I/O methods
  void write (ostream& out) const;
  void read (istream& in);
};

#endif /* MATRIX_ARRAY */
