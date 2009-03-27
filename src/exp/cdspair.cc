#include "cdspair.h"

CDS_pair_params::CDS_pair_params (PScores& pscore) :
  pscore_ptr (&pscore),
  first_pgroup (pscore.groups()),

  nullEmit (pscore.new_alphabet_group (CFG_alphabet, "nullEmit")),
  nullExtend (pscore.new_boolean_group ("nullExtend")),

  codEmit (pscore.new_alphabet_group (CFG_alphabet, 3, "codEmit")),
  codSub (CFG_alphabet_size * CFG_alphabet_size * CFG_alphabet_size),
  codGapOpen (pscore.new_boolean_group ("codGapOpen")),
  codGapExtend (pscore.new_boolean_group ("codGapExtend")),
  intron (pscore.new_boolean_group ("intron")),
  intronCons (pscore.new_boolean_group ("intronCons")),
  intExtend (pscore.new_boolean_group ("intExtend")),
  nucSub (CFG_alphabet_size),
  intGapOpen (pscore.new_boolean_group ("intGapOpen")),
  intGapExtend (pscore.new_boolean_group ("intGapExtend")),
  geneExtend (pscore.new_boolean_group ("geneExtend")),
  padOpen (pscore.new_boolean_group ("padOpen")),
  padExtend (pscore.new_boolean_group ("padExtend")),
  padGapOpen (pscore.new_boolean_group ("padGapOpen")),
  padGapExtend (pscore.new_boolean_group ("padGapExtend"))

{
  // initialize substitution PGroups
  for (int i = 0; i < (int) codSub.size(); ++i)
    {
      sstring name;
      name << "codSub" << codEmit.index2word(i);
      codSub[i] = pscore.new_alphabet_group (CFG_alphabet, 3, name.c_str());
    }

  for (int i = 0; i < (int) nucSub.size(); ++i)
    {
      sstring name;
      name << "nucSub" << CFG_alphabet.int2char(i);
      nucSub[i] = pscore.new_alphabet_group (CFG_alphabet, 1, name.c_str());
    }

  // set last_pgroup
  last_pgroup = pscore.groups() - 1;
}

Dirichlet_prior CDS_pair_params::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  // create the prior, assign idiotically simple pseudocounts to all PGroup's
  // (not quite Laplace, because k<1; but close)
  Dirichlet_prior prior (*pscore_ptr);
  for (int g = first_pgroup; g <= last_pgroup; ++g)
    prior.assign_Laplace (PGroup (g, pscore_ptr->group_size (g)), k);
  return prior;
}

int CDS_pair_state_enum::lex2 (int state, int pos, const vector<int>& x, const vector<int>& y)
{
  return
    pos == 0 ? state :
    (pos == 1 ? state + CONSSTATES * (1 + x[0] + 4 * y[0]) :
     (pos == 2 ? state + CONSSTATES * (17 + x[0] + 4 * (x[1] + 4 * (y[0] + 4 * y[1]))) :
      -1));
}

int CDS_pair_state_enum::lex (int state, int pos, int seq, const vector<int>& x)
{
  return
    pos == 0 ? state + TOTALCONSSTATES :
    (pos == 1 ? state + TOTALCONSSTATES + CODSTATES * (1 + seq + 2*x[0]) :
     (pos == 2 ? state + TOTALCONSSTATES + CODSTATES * (9 + seq + 2*x[0] + 8*x[1]) :
      -1));
}

CDS_pair::CDS_pair (PScores& pscore) :
  CDS_pair_params (pscore),
  Pair_PHMM (TotalStates, CFG_alphabet)
{
  // some useful variables
  const PFunc zero = 0.;
  const PFunc one = 1.;

  const int A = CFG_alphabet.char2int_strict ('A');
  // const int C = CFG_alphabet.char2int_strict ('C');
  const int G = CFG_alphabet.char2int_strict ('G');
  const int T = CFG_alphabet.char2int_strict ('T');

  const PFunc nullA = one / nullEmit[A];
  // const PFunc nullC = one / nullEmit[C];
  const PFunc nullG = one / nullEmit[G];
  const PFunc nullT = one / nullEmit[T];

  // initialize unlexed states
  // left-pad region
  init_emit (LPADMAT, EmitXY, zero);
  init_emit (LPADDEL, EmitX, one);
  init_emit (LPADINS, EmitY, one);
  state_type[LPADEND] = Null;

  transition (Start, LPADEND) = padOpen.n;
  transition (Start, LPADINS) = padOpen.y * padGapOpen.y;
  transition (Start, LPADDEL) = padOpen.y * padGapOpen.n * padGapOpen.y;
  transition (Start, LPADMAT) = padOpen.y * padGapOpen.n * padGapOpen.n;

  transition (LPADMAT, LPADINS) = padGapOpen.y;
  transition (LPADMAT, LPADDEL) = padGapOpen.n * padGapOpen.y;
  transition (LPADMAT, LPADMAT) = padGapOpen.n * padGapOpen.n * padExtend.y;
  transition (LPADMAT, LPADEND) = padGapOpen.n * padGapOpen.n * padExtend.n;

  transition (LPADINS, LPADINS) = padGapExtend.y;
  transition (LPADINS, LPADDEL) = padGapExtend.n * padGapOpen.y;
  transition (LPADINS, LPADMAT) = padGapExtend.n * padGapOpen.n * padExtend.y;
  transition (LPADINS, LPADEND) = padGapExtend.n * padGapOpen.n * padExtend.n;

  transition (LPADDEL, LPADDEL) = padGapExtend.y;
  transition (LPADDEL, LPADMAT) = padGapExtend.n * padExtend.y;
  transition (LPADDEL, LPADEND) = padGapExtend.n * padExtend.n;

  // right-pad region
  state_type[RPADSTART] = Null;
  init_emit (RPADMAT, EmitXY, zero);
  init_emit (RPADDEL, EmitX, one);
  init_emit (RPADINS, EmitY, one);

  transition (RPADSTART, End) = padOpen.n;
  transition (RPADSTART, RPADINS) = padOpen.y * padGapOpen.y;
  transition (RPADSTART, RPADDEL) = padOpen.y * padGapOpen.n * padGapOpen.y;
  transition (RPADSTART, RPADMAT) = padOpen.y * padGapOpen.n * padGapOpen.n;

  transition (RPADMAT, RPADINS) = padGapOpen.y;
  transition (RPADMAT, RPADDEL) = padGapOpen.n * padGapOpen.y;
  transition (RPADMAT, RPADMAT) = padGapOpen.n * padGapOpen.n * padExtend.y;
  transition (RPADMAT, End) = padGapOpen.n * padGapOpen.n * padExtend.n;

  transition (RPADINS, RPADINS) = padGapExtend.y;
  transition (RPADINS, RPADDEL) = padGapExtend.n * padGapOpen.y;
  transition (RPADINS, RPADMAT) = padGapExtend.n * padGapOpen.n * padExtend.y;
  transition (RPADINS, End) = padGapExtend.n * padGapOpen.n * padExtend.n;

  transition (RPADDEL, RPADDEL) = padGapExtend.y;
  transition (RPADDEL, RPADMAT) = padGapExtend.n * padExtend.y;
  transition (RPADDEL, End) = padGapExtend.n * padExtend.n;

  // emissions for pad regions
  for (int x = 0; x < CFG_alphabet_size; ++x)
    for (int y = 0; y < CFG_alphabet_size; ++y)
      pair_emit[LPADMAT](x,y) = pair_emit[RPADMAT](x,y) = nucSub[x][y] / nullEmit[y];

  // initialize lexed states
  vector<int> xcod (3);
  vector<int> ycod (3);
  for (xcod[0] = 0; xcod[0] < CFG_alphabet_size; ++xcod[0])
    for (xcod[1] = 0; xcod[1] < CFG_alphabet_size; ++xcod[1])
      for (ycod[0] = 0; ycod[0] < CFG_alphabet_size; ++ycod[0])
	for (ycod[1] = 0; ycod[1] < CFG_alphabet_size; ++ycod[1])
	  {
	    for (int pos = 0; pos < 3; ++pos)
	      {
		// conserved-codon state types
		const int i = lex2 (0, pos, xcod, ycod);

		init_emit (i+CONSSTART, pos<1 ? Null : EmitXY, zero);
		init_emit (i+CONSEND, pos<2 ? Null : EmitXY, zero);
		init_emit (i+CONSXYINT1, EmitXY, zero);
		init_emit (i+CONSXYINT2, EmitXY, zero);
		init_emit (i+CONSXYINTMAT, EmitXY, zero);
		init_emit (i+CONSXYINTDEL, EmitX, one);
		init_emit (i+CONSXYINTINS, EmitY, one);
		init_emit (i+CONSXYINT3, EmitXY, zero);
		init_emit (i+CONSXYINT4, EmitXY, zero);
		
		init_emit (i+CONSXINT1, EmitX, zero);
		init_emit (i+CONSXINT2, EmitX, zero);
		init_emit (i+CONSXINT, EmitX, one);
		init_emit (i+CONSXINT3, EmitX, zero);
		init_emit (i+CONSXINT4, EmitX, zero);

		init_emit (i+CONSYINT1, EmitY, zero);
		init_emit (i+CONSYINT2, EmitY, zero);
		init_emit (i+CONSYINT, EmitY, one);
		init_emit (i+CONSYINT3, EmitY, zero);
		init_emit (i+CONSYINT4, EmitY, zero);

		init_emit (i+CONSYINT1, EmitY, zero);
		init_emit (i+CONSYINT2, EmitY, zero);
		init_emit (i+CONSYINT, EmitY, one);
		init_emit (i+CONSYINT3, EmitY, zero);
		init_emit (i+CONSYINT4, EmitY, zero);

		// unconserved-codon state types
		for (int seq = 0; seq < 2; ++seq)
		  {
		    const int j = lex (0, pos, seq, xcod);
		    const State_type type = seq==0 ? EmitX : EmitY;

		    init_emit (j+CODSTART, pos<1 ? Null : type, zero);
		    init_emit (j+CODEND, pos<2 ? Null : type, zero);
		    init_emit (j+CODINT1, type, zero);
		    init_emit (j+CODINT2, type, zero);
		    init_emit (j+CODINT, type, one);
		    init_emit (j+CODINT3, type, zero);
		    init_emit (j+CODINT4, type, zero);
		  }
	      }

	    // intron emissions
	    for (int pos = 0; pos < 3; ++pos)
	      {
		// conserved-codon intron emissions
		const int i = lex2 (0, pos, xcod, ycod);

		pair_emit[i+CONSXYINT1](G,G) = nullG * nullG;
		pair_emit[i+CONSXYINT2](T,T) = nullT * nullT;
		pair_emit[i+CONSXYINT3](A,A) = nullA * nullA;
		pair_emit[i+CONSXYINT4](G,G) = nullG * nullG;

		for (int x = 0; x < CFG_alphabet_size; ++x)
		  for (int y = 0; y < CFG_alphabet_size; ++y)
		    pair_emit[i+CONSXYINTMAT](x,y) = nucSub[x][y] / nullEmit[y];

		single_emit[i+CONSXINT1][G] = nullG;
		single_emit[i+CONSXINT2][T] = nullT;
		single_emit[i+CONSXINT3][A] = nullA;
		single_emit[i+CONSXINT4][G] = nullG;

		single_emit[i+CONSYINT1][G] = nullG;
		single_emit[i+CONSYINT2][T] = nullT;
		single_emit[i+CONSYINT3][A] = nullA;
		single_emit[i+CONSYINT4][G] = nullG;

		// unconserved-codon intron emissions
		for (int seq = 0; seq < 2; ++seq)
		  {
		    const int j = lex (0, pos, seq, xcod);

		    single_emit[j+CODINT1][G] = nullG;
		    single_emit[j+CODINT2][T] = nullT;
		    single_emit[j+CODINT3][A] = nullA;
		    single_emit[j+CODINT4][G] = nullG;
		  }
	      }

	    // conserved codon emissions
	    for (xcod[2] = 0; xcod[2] < CFG_alphabet_size; ++xcod[2])
	      for (ycod[2] = 0; ycod[2] < CFG_alphabet_size; ++ycod[2])
		{
		  const int xc = codEmit.intvec2index (xcod);
		  const int yc = codEmit.intvec2index (ycod);

		  PFunc f = codEmit[xc] * codSub[xc][yc];
		  for (int k = 0; k < 3; ++k)
		    f /= nullEmit[xcod[k]] * nullEmit[ycod[k]];

		  pair_emit[lex2(CONSSTART,1,xcod,ycod)](xcod[0],ycod[0]) = one;
		  pair_emit[lex2(CONSSTART,2,xcod,ycod)](xcod[1],ycod[1]) = one;
		  pair_emit[lex2(CONSEND,2,xcod,ycod)](xcod[2],ycod[2]) = f;
		}

	    // unconserved codon emissions
	    for (int seq = 0; seq < 2; ++seq)
	      for (xcod[2] = 0; xcod[2] < CFG_alphabet_size; ++xcod[2])
		{
		  const int xc = codEmit.intvec2index (xcod);

		  PFunc f = codEmit[xc];
		  for (int k = 0; k < 3; ++k)
		    f /= nullEmit[xcod[k]];

		  single_emit[lex(CODSTART,1,seq,xcod)][xcod[0]] = one;
		  single_emit[lex(CODSTART,2,seq,xcod)][xcod[1]] = one;
		  single_emit[lex(CODEND,2,seq,xcod)][xcod[2]] = f;
		}

	    // transitions
	    for (int pos = 0; pos < 3; ++pos)
	      {
		// conserved codons
		const int i = lex2 (0, pos, xcod, ycod);

		// transitions within conserved codons
		if (pos > 1)
		  transition (lex2(CONSEND,pos-1,xcod,ycod), i+CONSSTART) = 1;

		// conserved-codon introns
		transition (i+CONSSTART, i+CONSEND) = intron.n;
		transition (i+CONSSTART, i+CONSXYINT1) = intron.y * intronCons.y;
		transition (i+CONSSTART, i+CONSXINT1) = intron.y * intronCons.n / 2;
		transition (i+CONSSTART, i+CONSYINT1) = intron.y * intronCons.n / 2;

		transition (i+CONSXYINT1, i+CONSXYINT2) = 1;
		transition (i+CONSXYINT2, i+CONSXYINTMAT) = 1;

		transition (i+CONSXYINTMAT, i+CONSXYINT3) = intGapOpen.n * intGapOpen.n * intExtend.n;
		transition (i+CONSXYINTMAT, i+CONSXYINTMAT) = intGapOpen.n * intGapOpen.n * intExtend.y;
		transition (i+CONSXYINTMAT, i+CONSXYINTINS) = intGapOpen.y;
		transition (i+CONSXYINTMAT, i+CONSXYINTDEL) = intGapOpen.n * intGapOpen.y;

		transition (i+CONSXYINTINS, i+CONSXYINT3) = intGapExtend.n * intGapOpen.n * intExtend.n;
		transition (i+CONSXYINTINS, i+CONSXYINTMAT) = intGapExtend.n * intGapOpen.n * intExtend.y;
		transition (i+CONSXYINTINS, i+CONSXYINTINS) = intGapExtend.y;
		transition (i+CONSXYINTINS, i+CONSXYINTDEL) = intGapExtend.n * intGapOpen.y;

		transition (i+CONSXYINTDEL, i+CONSXYINT3) = intGapExtend.n * intExtend.n;
		transition (i+CONSXYINTDEL, i+CONSXYINTMAT) = intGapExtend.n * intExtend.y;
		transition (i+CONSXYINTDEL, i+CONSXYINTDEL) = intGapExtend.y;

		transition (i+CONSXYINT3, i+CONSXYINT4) = 1;
		transition (i+CONSXYINT4, i+CONSEND) = 1;

		transition (i+CONSXINT1, i+CONSXINT2) = 1;
		transition (i+CONSXINT2, i+CONSXINT) = 1;
		transition (i+CONSXINT, i+CONSXINT) = intExtend.y;
		transition (i+CONSXINT, i+CONSXINT3) = intExtend.n;
		transition (i+CONSXINT3, i+CONSXINT4) = 1;
		transition (i+CONSXINT4, i+CONSEND) = 1;

		transition (i+CONSYINT1, i+CONSYINT2) = 1;
		transition (i+CONSYINT2, i+CONSYINT) = 1;
		transition (i+CONSYINT, i+CONSYINT) = intExtend.y;
		transition (i+CONSYINT, i+CONSYINT3) = intExtend.n;
		transition (i+CONSYINT3, i+CONSYINT4) = 1;
		transition (i+CONSYINT4, i+CONSEND) = 1;

		for (int seq = 0; seq < 2; ++seq)
		  {
		    // unconserved codons
		    const int j = lex (0, pos, seq, xcod);

		    // transitions within unconserved codons
		    if (pos > 1)
		      transition (lex(CODEND,pos-1,seq,xcod), j+CODSTART) = 1;

		    // unconserved-codon introns
		    transition (j+CODSTART, j+CODEND) = intron.n;
		    transition (j+CODSTART, j+CODINT1) = intron.y;
		    transition (j+CODINT1, j+CODINT2) = 1;
		    transition (j+CODINT2, j+CODINT) = 1;
		    transition (j+CODINT, j+CODINT) = intExtend.y;
		    transition (j+CODINT, j+CODINT3) = intExtend.n;
		    transition (j+CODINT3, j+CODINT4) = 1;
		    transition (j+CODINT4, j+CODEND) = 1;
		  }
	      }

	    // transitions between codons and LPADEND/RPADSTART
	    transition (LPADEND, lex(CODSTART,0,0,xcod)) = codGapOpen.y;
	    transition (LPADEND, lex(CODSTART,0,1,ycod)) = codGapOpen.n * codGapOpen.y;
	    transition (LPADEND, lex2(CODSTART,0,xcod,ycod)) = codGapOpen.n * codGapOpen.n * geneExtend.y;
	    transition (LPADEND, RPADSTART) = codGapOpen.n * codGapOpen.n * geneExtend.n;

	    transition (lex2(CODEND,2,xcod,ycod), lex(CODSTART,0,0,xcod)) = codGapOpen.y;
	    transition (lex2(CODEND,2,xcod,ycod), lex(CODSTART,0,1,ycod)) = codGapOpen.n * codGapOpen.y;
	    transition (lex2(CODEND,2,xcod,ycod), lex2(CODSTART,0,xcod,ycod)) = codGapOpen.n * codGapOpen.n * geneExtend.y;
	    transition (lex2(CODEND,2,xcod,ycod), RPADSTART) = codGapOpen.n * codGapOpen.n * geneExtend.n;

	    transition (lex(CODEND,2,0,xcod), lex(CODSTART,0,0,xcod)) = codGapExtend.y;
	    transition (lex(CODEND,2,0,xcod), lex(CODSTART,0,1,ycod)) = codGapExtend.n * codGapOpen.y;
	    transition (lex(CODEND,2,0,xcod), lex2(CODSTART,0,xcod,ycod)) = codGapExtend.n * codGapOpen.n * geneExtend.y;
	    transition (lex(CODEND,2,0,xcod), RPADSTART) = codGapExtend.n * codGapOpen.n * geneExtend.n;

	    transition (lex(CODEND,2,1,ycod), lex(CODSTART,0,1,ycod)) = codGapExtend.y;
	    transition (lex(CODEND,2,1,ycod), lex2(CODSTART,0,xcod,ycod)) = codGapExtend.n * geneExtend.y;
	    transition (lex(CODEND,2,1,ycod), RPADSTART) = codGapExtend.n * geneExtend.n;
	  }

  // factor in null_extend
  const PFunc g = nullExtend.y * nullExtend.y;
  const PFunc h = nullExtend.n * nullExtend.n;
  for (int s = 0; s < states(); ++s)
    {
      if (state_type[s] == EmitXY)
	for_contents (array2d<PFunc>, pair_emit[s], f)
	  if (!f->is_zero())
	    *f /= g;
      else if (state_type[s] != Null)
	for_contents (vector<PFunc>, single_emit[s], f)
	  if (!f->is_zero())
	    *f /= nullExtend.y;
      PFunc& f = transition (s, End);
      if (!f.is_zero())
	f /= h;
    }
}
