#ifndef SHORTLIST_INCLUDED
#define SHORTLIST_INCLUDED

#include "util/array2d.h"
#include "empath/trainer.h"
#include "empath/multihit.h"

// Seq_pos: sequence index + position
struct Seq_pos
{
  int seq;
  int pos;
  Seq_pos (int seq, int pos) : seq (seq), pos (pos) { }
  inline bool operator== (const Seq_pos& sp) const { return seq == sp.seq && pos == sp.pos; }
  inline bool operator< (const Seq_pos& sp) const { return seq < sp.seq || (seq == sp.seq && pos < sp.pos); }
};

// sort function for Seq_pos's
struct Seq_pos_by_likelihood : binary_function <const Seq_pos&, const Seq_pos&, bool>
{
  const vector<vector<Loge> >& log_like;
  Seq_pos_by_likelihood (const vector<vector<Loge> >& log_like) : log_like (log_like) { }
  inline bool operator() (const Seq_pos& a, const Seq_pos& b) const
  { return log_like[a.seq][a.pos] > log_like[b.seq][b.pos]; }
};

// Shortlist: list of candidate Seq_pos seeds
struct Shortlist : Piper
{
  // member data
  sstring name;  // name of this shortlist
  const FASTA_sequence_database& seq_db;  // the database
  const Sequence_database_index& index;   // the database index
  vector<Seq_pos> seed;  // ordered seed list (best first)
  vector<vector<Loge> > loglike;  // a log-likelihood for every seed
  Seq_pos_by_likelihood by_loglike;  // the sort function
  const bool diagnostic_flag;  // flag to run diagnostic checks
  // constructor, virtual destructor
  Shortlist (const FASTA_sequence_database& seq_db, int max_fork);
  virtual ~Shortlist();
  // methods to allocate & free loglike table
  void alloc_loglike();
  void free_loglike();
  // seeding methods
  void seed_everywhere (int motif_len);  // seeds everywhere
  void seed_locally_repeated (int motif_len, int word_len, int min_reps);  // seeds on locally repeated words
  void seed_globally_repeated (int motif_len, int word_len, int min_reps);  // seeds on globally repeated words
  void copy_seeds (const Shortlist& source, int N);  // copies first N seeds from source
  // full sort method: rebuilds entire loglike table
  void full_sort();
  // partial resort method: re-evaluates top N seeds, re-sorts, iterates until top N are stable
  // ...or until all with log-likelihood >= min_loglike are stable (in which case, set N=0)
  void partial_sort (int N = 0, Loge min_loglike = -InfinityLoge);
  // method to return the number of seeds above a given log-likelihood, following a sort
  int n_above (Loge min_loglike) const;
  // virtual method to calculate log-likelihood for any Seq_pos
  virtual Loge compute_loglike (const Seq_pos& seqpos) = 0;
  // virtual method to reset any internal cache, called by full_sort() and partial_sort()
  virtual void reset_cache();
  // virtual debugging/diagnostic method, called if a seed likelihood increases, or if "-log SHORTLIST_DIAGNOSTIC" specified
  virtual void run_diagnostic (const Seq_pos& seqpos, Loge loglike, bool be_verbose);
  // output method
  void show (ostream& o) const;
};

// Seed_drill: chain of Shortlists along which seeds are passed
// courtesy of Jethro Tull
class Seed_drill
{
private:
  bool first_pass;  // TRUE until first call to refresh()
  vector<Shortlist*> shortlist;  // shortlists
  vector<int> len;  // number of seeds to keep from each shortlist
  vector<Loge> min_loglike;  // minimum log-likelihood of seeds to keep from each shortlist
public:
  Seed_drill();
  void clear();
  void delete_shortlists();  // calls clear()
  void add_fixed_length_shortlist (Shortlist* s, int s_len);
  void add_variable_length_shortlist (Shortlist* s, Loge s_loglike);
  Shortlist& front() { return *shortlist.front(); }
  void refresh();  // sorts each shortlist in turn, passing seeds along
  const Seq_pos& best_seed() const;
};

struct Self_similarity_shortlist : Shortlist
{
  // cache
  vector<vector<int> > first_instance;  // table of indices to first instance of each seed word within a sequence
  map<Seq_pos,Score> fresh;  // set of freshly computed seed likelihoods
  // logging flags
  bool dp_logging;
  bool long_dp_logging;
  // data
  const Alphabet& alphabet;
  const int motif_len;
  const int offset;
  const int mask_metascore_idx;
  const bool revcomp;
  array2d<Score> submat;
  vector<Score> max_sym_sc;
  Digitized_biosequence revcomp_target;  // temporary variable for revcomp-ing
  // methods
  void reset_cache();  // clears fresh
  Loge compute_loglike (const Seq_pos& seqpos);
  Self_similarity_shortlist (const FASTA_sequence_database& seq_db,
			     int max_fork,
			     const Alphabet& alphabet,
			     int motif_len,
			     int offset,
			     Prob p_match,
			     bool revcomp,
			     int mask_metascore_idx);
  // slow debugging/diagnostic methods for computing log-likelihood
  void run_diagnostic (const Seq_pos& seqpos, Loge loglike, bool be_verbose);  // compares loglike with compute_loglike_slow()
  Loge compute_loglike_slow (const Seq_pos& seqpos);  // doesn't attempt to accelerate DP
  Score compute_match_sc (const int* target, const int* query, const Score* query_mask) const;
  void show_dp_vec (const vector<Score>& fwd, const vector<Score>& rev, ostream& out) const;
};

struct Local_shortlist : Shortlist
{
  // data
  Local_trainer& trainer;
  // methods
  Loge compute_loglike (const Seq_pos& seqpos);
  Local_shortlist (const FASTA_sequence_database& seq_db, int max_fork, Local_trainer& trainer);
};

struct Multihit_shortlist : Shortlist
{
  // data
  Multihitter& model;
  // the following vector gives, for each seq, a list of similar seqs according to distmat
  // if this vector is empty, then the full database is screened on every compute_loglike()
  vector<vector<Named_profile*> > similar_np;
  // methods
  Loge compute_loglike (const Seq_pos& seqpos);
  // constructors
  Multihit_shortlist (const FASTA_sequence_database& seq_db, Multihitter& model);
  Multihit_shortlist (const FASTA_sequence_database& seq_db, Multihitter& model, Dist_func& dist, int n_similar);
};

#endif /* SHORTLIST_INCLUDED */
