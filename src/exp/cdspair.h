#ifndef CDSPAIR_INCLUDED
#define CDSPAIR_INCLUDED

#include "hmm/pairphmm.h"
#include "scfg/paircfg.h"
#include "seq/dirichlet.h"

class CDS_pair_params
{
private:
  // data
  PScores* pscore_ptr;
  int first_pgroup;   // index of first PGroup created by this object
  int last_pgroup;    // index of last PGroup

public:
  // Probability parameters
  Alphabet_group nullEmit;        // null model (also intron model)
  Boolean_group nullExtend;       // null model extend probability
  Alphabet_group codEmit;         // equilibrium distribution over codons
  vector<Alphabet_group> codSub;  // conditional codon substitution probabilities
  Boolean_group codGapOpen;       // gap open probability in coding sequence
  Boolean_group codGapExtend;     // gap extend probability in coding sequence
  Boolean_group intron;           // intron probability
  Boolean_group intronCons;       // probability that an intron is conserved
  Boolean_group intExtend;        // intron extend probability
  vector<Alphabet_group> nucSub;  // conditional nucleotide substitution probabilities (introns & padding regions)
  Boolean_group intGapOpen;       // gap open probability in introns
  Boolean_group intGapExtend;     // gap extend probability in introns
  Boolean_group geneExtend;       // gene extension probability
  Boolean_group padOpen;          // padding region open probability
  Boolean_group padExtend;        // padding region extension probability
  Boolean_group padGapOpen;       // padding region gap-open probability
  Boolean_group padGapExtend;     // padding region gap-extension probability

  // constructor
  CDS_pair_params (PScores& pscore);

  // method to generate prior
  Dirichlet_prior default_prior() const;
};

struct CDS_pair_state_enum
{
  enum {
    // conserved codon states (1 0-lexed, 16 2-lexed and 256 4-lexed)
    CONSSTART = 0,// type EmitXY (positions 2,3) or Null (position 1), conserved codon position
    CONSEND,      // type EmitXY (position 3) or Null (positions 1,2), conserved codon position
    CONSXYINT1,   // type EmitXY, conserved intron, emits first base (G)
    CONSXYINT2,   // type EmitXY, conserved intron, emits second base (T)
    CONSXYINTMAT, // type EmitXY, conserved intron, match state
    CONSXYINTDEL, // type EmitX, conserved intron, delete state
    CONSXYINTINS, // type EmitY, conserved intron, insert state
    CONSXYINT3,   // type EmitXY, conserved intron, penultimate base (A)
    CONSXYINT4,   // type EmitXY, conserved intron, last base (G)

    CONSXINT1,    // type EmitX, deleted intron, emits first base (G)
    CONSXINT2,    // type EmitX, deleted intron, emits second base (T)
    CONSXINT,     // type EmitX, deleted intron, emit state
    CONSXINT3,    // type EmitX, deleted intron, penultimate base (A)
    CONSXINT4,    // type EmitX, deleted intron, last base (G)

    CONSYINT1,    // type EmitY, inserted intron, emits first base (G)
    CONSYINT2,    // type EmitY, inserted intron, emits second base (T)
    CONSYINT,     // type EmitY, inserted intron, emit state
    CONSYINT3,    // type EmitY, inserted intron, penultimate base (A)
    CONSYINT4,    // type EmitY, inserted intron, last base (G)

    CONSSTATES = 19,        // number of conserved codon-related states per lex-context
    TOTALCONSSTATES = 5187, // total number of conserved codon-related states = 19*(1+16+256)

    // unconserved codon states (1 0-lexed, 4 1-lexed and 16 2-lexed)
    // for each, there is one EmitX & one EmitY
    CODSTART = 0, // type Emit (positions 2,3) or Null (position 1), unconserved codon position
    CODEND,       // type Emit (position 3) or Null (positions 1,2), unconserved codon position
    CODINT1,      // unconserved intron, emits first base (G)
    CODINT2,      // unconserved intron, emits second base (T)
    CODINT,       // unconserved intron, emit state
    CODINT3,      // unconserved intron, penultimate base (A)
    CODINT4,      // unconserved intron, last base (G)

    CODSTATES = 7,        // number of conserved codon-related states per lex-context per sequence
    TOTALCODSTATES = 294, // total number of conserved codon-related states = 7*2*(1+4+16)

    // Unlexed states
    // TOTALCONSSTATES + TOTALCODSTATES = 5481
    LPADMAT = 5481,  // type EmitXY, left-pad match state
    LPADDEL,         // type EmitX, left-pad delete state
    LPADINS,         // type EmitY, left-pad insert state
    LPADEND,         // type Null, end of left-pad region, pseudo-Start state
    RPADSTART,       // type Null, start of right-pad region, pseudo-End state
    RPADMAT,         // type EmitXY, right-pad match state
    RPADDEL,         // type EmitX, right-pad delete state
    RPADINS,         // type EmitY, right-pad insert state

    TotalStates      // yields total number of states
  };

  // lexed state index helpers
  static int lex2 (int state, int pos, const vector<int>& x, const vector<int>& y);  // for (2*pos)-lexed states
  static int lex (int state, int pos, int seq, const vector<int>& x);  // for (pos)-lexed states, with seq = 0(X) or 1(Y)
};

struct CDS_pair : CDS_pair_state_enum, CDS_pair_params, Pair_PHMM
{
  // constructor
  CDS_pair (PScores& pscore);
};

#endif /* CDSPAIR_INCLUDED */
