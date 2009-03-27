#ifndef DIRICHLET_INCLUDED
#define DIRICHLET_INCLUDED

#include <list>
#include <algorithm>
#include "seq/pvar.h"

// Multigroup: a list of PGroups
// some day, may integrate this into Dirichlet_assignment etc
struct Multigroup : vector<PGroup>
{
  Multigroup() : vector<PGroup>() { }
  Multigroup (PGroup g) : vector<PGroup>() { push_back(g); }
  Multigroup (const vector<PGroup>& v) : vector<PGroup> (v) { }

  Multigroup& operator*= (PGroup g)     { push_back(g); return *this; }
  Multigroup& operator*= (Multigroup m) { for_const_contents (Multigroup, m, g) push_back(*g); return *this; }

  friend Multigroup operator* (PGroup l,     PGroup r)     { Multigroup m; m *= l; m *= r; return m; }
  friend Multigroup operator* (PGroup l,     Multigroup r) { r.insert (r.begin(), l); return r; }
  friend Multigroup operator* (Multigroup l, PGroup r)     { l *= r; return l; }
  friend Multigroup operator* (Multigroup l, Multigroup r) { l *= r; return l; }
};

Multigroup operator* (PGroup l,     PGroup r);
Multigroup operator* (PGroup l,     Multigroup r);
Multigroup operator* (Multigroup l, PGroup r);
Multigroup operator* (Multigroup l, Multigroup r); 

// Dirichlet_mixture class represents a mixture prior for a Multigroup
struct Dirichlet_mixture
{
  vector<Loge>                   log_cpt_prior;
  vector<vector<vector<Prob> > > alpha;   // alpha[component][group_idx][var_idx]
  
  int components() const { return alpha.size(); }
  int groups() const { return alpha[0].size(); }
  int vars (int group_idx) const { return alpha[0][group_idx].size(); }

  void assert_consistent() const;   // checks that alpha[][][] is rectangular

  // of the following three constructors, only the first initialises the alphas to anything except zero:

  Dirichlet_mixture (const vector<Prob>& pseudocounts);               // single-component single-group constructor
  Dirichlet_mixture (int group_size, int components);                 // multi-component single-group constructor
  Dirichlet_mixture (const vector<int>& group_size, int components);  // multi-component multi-group constructor
  Dirichlet_mixture (const Multigroup& multigroup, int components);   // multi-component multi-group constructor

  // methods for doing actual math
  // always call calc_log_evidence BEFORE postcalc_optimise_scores

  void calc_log_evidence (const PCounts& pcounts,
			  const vector<PGroup>& pgroups,
			  vector<Loge>& log_evidence,
			  Loge& log_total_evidence) const;

  void postcalc_optimise_scores (const PCounts& pcounts,
				 const vector<PGroup>& pgroups,
				 const vector<Loge>& log_evidence,
				 const Loge& log_total_evidence,
				 PScores& pscores) const;

  void initialise_default_scores (PScores& pscores,
				  const vector<PGroup>& pgroups) const;
};

// Laplace prior: +1 pseudocounts (actually, can be +k, where k is any real)
// The krange parameter allows the pseudocounts to increase slightly (if k>1) or taper off (if k<1) over the PGroup.
// The pseudocount for PVar V in a PGroup of size G (with 0<=V<G) is k*pow(krange,V/(G-1))
struct Laplace_prior : Dirichlet_mixture
{
  Laplace_prior (const Multigroup& multigroup, double k = 1., double krange = 1.);
  Laplace_prior (const PGroup& pgroup, double k = 1., double krange = 1.);
};

// Alphabet prior, for a Multigroup consisting of a single Alphabet_group
struct Alphabet_prior : Dirichlet_mixture
{
  const Alphabet& alphabet;
  Alphabet_prior (const Alphabet& alphabet, int components);
};

// Kimmen Sjolander's 9-component protein priors
struct Sjolander_prior : Alphabet_prior
{
  Sjolander_prior (const Alphabet& alphabet = Protein_alphabet);
};

struct Dirichlet_assignment
{
  // variables that specify this assignment

  Dirichlet_mixture mixture;
  Multigroup        pgroups;

  // temporary variables for evidence calculations

  vector<Loge> log_evidence;
  Loge         log_total_evidence;

  // constructors

  Dirichlet_assignment() : mixture(0,0), pgroups() { }
  Dirichlet_assignment (const Dirichlet_mixture& mixture);
  Dirichlet_assignment (const Dirichlet_mixture& mixture, const vector<PGroup>& pgroups);
  Dirichlet_assignment (const Dirichlet_assignment& a) : mixture (a.mixture), pgroups (a.pgroups) { }
  
  // copy operator
  Dirichlet_assignment& operator= (const Dirichlet_assignment& a) { mixture = a.mixture; pgroups = a.pgroups; return *this; }

  // bridge methods to Dirichlet_mixture
  // optimise_scores calls calc_log_evidence, so there is no need to call calc_log_evidence unless that's all you want to do

  void calc_log_evidence (const PCounts& pcounts);
  void optimise_scores (const PCounts& pcounts, PScores& pscores);
  void initialise_default_scores (PScores& pscores) const;
};

// Dirichlet prior
struct Dirichlet_prior
{
  // Strict-weak-ordering functor for testing whether two Dirichlet_assignments refer to the same Multigroup
  //  struct Compare_multigroups
  //  {
  //    bool operator() (const Multigroup& l, const Multigroup& r) const
  //    { return lexicographical_compare (l.begin(), l.end(), r.begin(), r.end()); }
  //  };

  // typedef for a set of Dirichlet_assignments
  typedef map <Multigroup, Dirichlet_assignment> Assignment_set;

  // data
  PScores* pscores;  // relying on this to point to a valid object is EXTREMELY risky. Best to set it yourself before calling optimise() or initialise()
  Assignment_set assignment;

  // constructors
  Dirichlet_prior() : pscores(0) { }
  Dirichlet_prior (PScores& pscores) : pscores (&pscores) { }

  // build methods
  void clear();
  void assign (const vector<PGroup>& pgroups, const Dirichlet_mixture& mixture);
  void assign (const PGroup& pgroup, const Dirichlet_mixture& mixture);
  void assign_Laplace (const PGroup& pgroup, const double laplace_pseudocount = 1., const double pseudocount_range_ratio = 1.);
  void assign_Laplace (PScores& pscores, const double laplace_pseudocount = 1., const double pseudocount_range_ratio = 1.);  // calls assign_Laplace for each PGroup
  void assign_Laplace (PScores& pscores, set<int>& mutable_pgroups, const double laplace_pseudocount = 1., const double pseudocount_range_ratio = 1.);  // calls assign_Laplace for each PGroup in mutable_pgroups

  void add (const Dirichlet_prior& prior);

  // application
  void optimise (const PCounts& pcounts);
  void initialise() const;
};

#endif
