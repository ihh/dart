#ifndef PROGHMM_INCLUDED
#define PROGHMM_INCLUDED

#include "hmm/pairphmm.h"
#include "seq/dirichlet.h"

struct Prog_HMM_base : Pair_PHMM
{
  // data
  PScores pscore;
  Dirichlet_prior prior;
  // constructor
  Prog_HMM_base (int states, const Alphabet& alphabet);
  // virtual destructor
  virtual ~Prog_HMM_base();
  // helpers
  void assign_Laplace_prior();
};

struct Prog_HMM : Prog_HMM_base
{
  // state enumeration
  enum { ZoneStartBase = 0, ZoneEndBase = 1, MatchBase = 2, InsBase = 3, DelBase = 4 };
  // data
  // number of indel types & zones
  int indelTypes, zones;
  // PGroup's
  vector<Boolean_group> stop_pg;  // indexed by zone
  vector<Boolean_group> indel_begin_pg;  // indexed by zone
  vector<PGroup> indel_type_pg;  // indexed by zone
  vector<vector<Boolean_group> > indel_extend_pg;  // indexed by zone, then indel type
  vector<vector<Boolean_group> > indel_swap_pg;  // indexed by zone, then indel type
  vector<Alphabet_group> single_pg;  // indexed by zone
  vector<vector<Alphabet_group> > sub_pg;  // indexed by zone, then source symbol
  PGroup start_zone_pg;
  vector<PGroup> next_zone_pg;  // indexed by source zone
  // constructor
  Prog_HMM (int indelTypes, int zones, const Alphabet& alphabet);
  // state accessors
  inline int ZoneStart (int zone) const;
  inline int ZoneEnd (int zone) const;
  inline int Match (int zone) const;
  inline int Ins (int type, int zone) const;
  inline int Del (int type, int zone) const;
  static inline int statesPerIndelType();
  static inline int statesPerZone (int indelTypes);
  static inline int calcTotalStates (int indelTypes, int zones);
};


// inline method defs

// Prog_HMM

int Prog_HMM::ZoneStart (int zone) const { return zone * statesPerZone(indelTypes) + ZoneStartBase; }
int Prog_HMM::ZoneEnd (int zone) const { return zone * statesPerZone(indelTypes) + ZoneEndBase; }
int Prog_HMM::Match (int zone) const { return zone * statesPerZone(indelTypes) + MatchBase; }
int Prog_HMM::Ins (int zone, int type) const { return zone * statesPerZone(indelTypes) + 2*type + InsBase; }
int Prog_HMM::Del (int zone, int type) const { return zone * statesPerZone(indelTypes) + 2*type + DelBase; }

int Prog_HMM::statesPerZone (int indelTypes) { return 2 * indelTypes + 3; }
int Prog_HMM::calcTotalStates (int indelTypes, int zones) { return zones * statesPerZone(indelTypes); }

#endif /* PROGHMM_INCLUDED */
