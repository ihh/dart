#ifndef PSWMODEL
#define PSWMODEL

#include "hmm/pairphmm.h"

class PSW_vars
{
public:
  PScope* scope;
  const Alphabet* alphabet;
  int time_bins;
  // Transition matrix for branch transducer
  //       I   W1    M   W2    D    E
  //  S    a !a!b    0  !ab    0    0
  //  I    c !c!d    0  !cd    0    0
  // W1    0    0    1    0    0    1
  //  M    e !e!f    0  !ef    0    0
  // W2    0    0    0    0    1    1
  //  D    g !g!h    0  !gh    0    0
  vector<Boolean_group> a, b, c, d, e, f, g, h;
  Boolean_group aroot, croot;  // infinite-time a & c
  Alphabet_group emit;
  vector<vector<Alphabet_group> > subst;  // indexed by time, then ancestral symbol
  // constructors, initialiser
  PSW_vars();
  PSW_vars (PScope& scope, const Alphabet& alphabet, int time_bins);
  void init (PScope& scope, const Alphabet& alphabet, int time_bins);
private:
  Boolean_group new_boolean_group (const char* prefix, int time);
  Alphabet_group new_alphabet_group (const char* prefix, int time);
};

struct Zoned_PSW_vars
{
  PScope& scope;
  const Alphabet& alphabet;
  const int zones;
  const int time_bins;
  vector<PSW_vars> zone_vars;  // indexed by zone
  PGroup zone_start;  // PGroup indexed by zone
  vector<PGroup> zone_trans;  // indexed by src zone; PGroup indexed by dest zone
  // constructor
  Zoned_PSW_vars (PScope& scope, const Alphabet& alphabet, int zones, int time_bins);
};

struct PSW_HMM_state_indices
{
  static inline int SS (int zone) { return zone*6 + 0; }   // state type: Null
  static inline int EE (int zone) { return zone*6 + 1; }   // state type: Null
  static inline int SI (int zone) { return zone*6 + 2; }   // state type: EmitY
  static inline int II (int zone) { return zone*6 + 3; }   // state type: EmitY
  static inline int IM (int zone) { return zone*6 + 4; }   // state type: EmitXY
  static inline int ID (int zone) { return zone*6 + 5; }   // state type: EmitX
  static inline int total_states (int zones) { return zones*6; }
};

struct PSW_PHMM : Pair_PHMM, PSW_HMM_state_indices
{
  const Zoned_PSW_vars& vars;
  const int time;
  const int zones;
  PSW_PHMM (const Zoned_PSW_vars& vars, int time, bool allow_null_cycles = FALSE);
};

#endif /* PSWMODEL */
