#ifndef TKFSEQUENCE_INCLUDED
#define TKFSEQUENCE_INCLUDED

#include <fstream>
#include <deque>
#include "tkf/tkfparams.h"
#include "util/score.h"

struct TKF_sequence : Digitized_biosequence
{
  Pairwise_path path;
  double        log_path_likelihood;
  double        log_emit_likelihood;

  vector<bool>& parent_path() { return path[0]; }
  vector<bool>& my_path()     { return path[1]; }

  void generate (const TKF_params& params);
  void generate_predicated_on_length (const TKF_params& params, int length);   // NB likelihood of expected length will still be factored into log_path_likelihood, so it's not conditioned (just predicated?)
  void generate_expected_length (const TKF_params& params) { generate_predicated_on_length (params, (int) (0.5 + params.expected_sequence_length())); }   // rounding up the length is gooey but works
  
  void evolve (const TKF_params& params, const Digitized_biosequence& parent, double time);
};



struct TKF_emit_align : Tree_alignment
{
  const TKF_params&         params;
  vector<TKF_sequence>      tkf;
  vector<Named_profile>     np;

  TKF_emit_align (const TKF_params& params) : params(params) { }

  void   emit (bool use_expected_len = 0);

  double log_path_likelihood() const;
  double log_emit_likelihood() const;
  double log_likelihood() const { return log_path_likelihood() + log_emit_likelihood(); }
  
  int score() const;
};



#endif
