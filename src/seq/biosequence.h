// biosequence code using STL
// ihh@lanl.gov, 1999/04/28

#ifndef BIOSEQUENCE_INCLUDED
#define BIOSEQUENCE_INCLUDED

#include <stdio.h>
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include "util/sstring.h"
#include "util/macros.h"
#include "util/score.h"
#include "util/phonebook.h"
#include "util/rnd.h"
#include "util/map_keys.h"

#define Alphabet_nonsense_char '?'

// a Biosequence is a text representation of a sequence
//
typedef sstring Biosequence;

// a Digitized_biosequence is a numerical representation of a sequence, related to a Biosequence by an Alphabet (see below)
// NB implementation isn't quite there (will throw an exception for negative (=ambiguous) symbols)
//
class Digitized_biosequence : public vector<int>
{
public:
  Digitized_biosequence() : vector<int>() { }
  Digitized_biosequence (const vector<int>& intvec) : vector<int> (intvec) { }
  Digitized_biosequence (size_type sz, int i = -1) : vector<int> (sz, i) { }
  int  null_score (const vector<int>& null_model) const;
  int  stuttered_null_score (const vector<int>& null_model, const int stutter_sc) const;
  void count_symbol_frequencies (vector<double>& composition) const;
  // the "context" methods a Digitized_biosequence where every site
  // encodes the last 'context_order' residues of the parent sequence
  // ...a hacky way of doing higher-order Markov chains. probably best to ignore these methods.
  Digitized_biosequence make_context_dsq (int alphabet_size, int context_order) const;
  Digitized_biosequence make_reverse_context_dsq (int alphabet_size, int context_order) const;
  static Digitized_biosequence context_to_nmer (int context, int alphabet_size, int context_order);
};

// a Weight_profile is an ungapped profile in probability space
// it should probably have all the same methods as Score_profile (see below), but it is more rarely used
//
typedef pair <int, Prob>                   Symbol_weight;
typedef map <int, Prob>                    Symbol_weight_map;
typedef Symbol_weight_map::const_iterator  Symbol_weight_iterator;

class Weight_profile : public vector <Symbol_weight_map>
{
 public:
  Weight_profile () {}
  Weight_profile (const Digitized_biosequence& dsq);
  Weight_profile (int len, const Symbol_weight_map& weight_map) : vector<Symbol_weight_map> (len, weight_map) {}

  static int weight_map_consensus (const Symbol_weight_map& w);           // returns most heavily weighted symbol in a Symbol_score_map
  int        consensus (int pos) const { return weight_map_consensus ((*this)[pos]); }

  Digitized_biosequence consensus_dsq() const;

  void       normalise();               // for each residue, re-scales symbol weights so that probability sum is 1
  void       count_symbol_frequencies (vector<double>& composition) const;
};

// a Vector_weight_profile is an ungapped probability profile like a Weight_profile,
// but stored as a vector
typedef vector<vector<Prob> > Vector_weight_profile;

// a Score_profile is an ungapped profile in score space
//
typedef pair <int, Score>                 Symbol_score;
typedef map <int, Score>                  Symbol_score_map;
typedef Symbol_score_map::const_iterator  Symbol_score_iterator;

class Score_profile : public vector <Symbol_score_map>
{
 public:
  Score_profile () {}
  Score_profile (const Digitized_biosequence& dsq);
  Score_profile (const Weight_profile& profile);
  Score_profile (int len, const Symbol_score_map& score_map) : vector<Symbol_score_map> (len, score_map) {}

  // Should really make Symbol_score_map a proper class and put this class's static methods there (ditto Symbol_weight_map)
  //
  static int     score_map_consensus (const Symbol_score_map& s);                // returns highest scoring symbol in a Symbol_score_map
  int            consensus (int pos) const { return score_map_consensus ((*this)[pos]); }

  static int     score_map_sample (const Symbol_score_map& s, double kT = 1);    // samples from a Symbol_score_map
  int            sample (int pos, double kT = 1) const { return score_map_sample ((*this)[pos]); }

  static int     score_map_score (const Symbol_score_map& ssm, const vector<int>& null_model);
  int            null_score (const vector<int>& null_model) const;

  void           normalise();               // for each residue, re-scales symbol scores so that probability sum is 1
  
  Digitized_biosequence consensus_dsq() const;
  Digitized_biosequence sample_dsq (double kT = 1) const;

  Score_profile  evolve (const array2d<int>& substitution_matrix) const;            // evolves the sequence (substitutions only)
  vector<int>    null_scores_by_residue (const vector<int>& null_model) const;      // returns null score for each residue in the sequence

  int            inner_product (const Score_profile& prof) const;    // score of product_{columns} sum_{residues in column} P(residue|this) . P(residue|prof)
  
  Weight_profile score2weight() const;         // turns this Score_profile into a Weight_profile

  void           count_symbol_frequencies (vector<double>& composition) const;       // count symbol frequencies
};

// Alphabet class
//
class Alphabet
{
  // every symbol in the alphabet corresponds to a non-negative integer 0,1,2...(size-1)
  // negative integers represent ambiguous symbols, e.g. "N" "Y" "R" for DNA, "X" for proteins, and wildcards "*"
  // these ambiguities can be represented by a weighted sum over the real symbols in the alphabet
  //  e.g. "Y" = 0.5 * "C" + 0.5 * "T" (NB this is a probability distribution)
  //       "*" = "A" + "C" + "G" + "T" (NB this is not a probability distribution)
  // such weighted sums are represented by Symbol_weight_map's (in probability space) or Symbol_score_map's (in score space)
  //  these representations can be resolved by calling int2weight() or int2score()
  //  in fact there are a million trillion zillion different conversion functions we have to define, and some are as yet unused & unimplemented.
  //
 protected:
  int alphabet_size;  // number of nondegenerate symbols

  // arrays indexed by the character representation of an alphabet token
  vector<int> char2int_table;     // a value of (alphabet_size) indicates unknown
  vector<double> char_aff_table;  // table of char affinities for this alphabet; used to match sequences to alphabets

  // arrays indexed by the integer representation of an alphabet token
  vector<int> nondegen_comp;      // complements for nondegenerate chars
  vector<int> degen_comp;         // complements for degenerate chars

  vector<sstring> tok;  // tokens for nondegenerate symbols (a token can be multiple chars)
  sstring chars_lc;  // nondegenerate chars, lower case
  sstring chars_uc;  // nondegenerate chars, upper case
  sstring degen_lc;  // degenerate chars, lower case
  sstring display_lc;  // what to output for nondegenerate symbols with no chars (used by hidden alphabets)

  vector<Symbol_weight_map> nondegen_weight;  // nondegenerate weight mappings
  vector<Symbol_score_map>  nondegen_score;   // nondegenerate score mappings
  vector<Symbol_weight_map> degen_weight;     // degenerate weight mappings
  vector<Symbol_score_map>  degen_score;      // degenerate score mappings

  int unknown_int;  // "N" for DNA, "X" for protein

  bool warn_randomise;  // flag saying whether to warn when randomising degenerate chars; set to 0 after first warning
  bool warn_unknown;    // flag for unknown chars warning

 public:
  // name
  sstring name;  // name of alphabet, e.g. "DNA", "Protein"

  // case-sensitivity flag
  bool case_sensitive;

  // multiplier for counts, to turn fractional degenerate counts into integers. 12 for DNA/RNA, 20 for protein
  //   (explanation: ambiguated counts for DNA can be 1/2 (e.g. "r", "y"), 1/3 (e.g. "h", "v"), 1/4 ("n")
  //    lowest common denominator of (1/2, 1/3, 1/4) is 12)
  unsigned int degenerate_lcd;

  // Felsenstein wildcards
  Symbol_score_map wild_ssm;
  Symbol_weight_map wild_swm;

  // builder methods
  void reset (const char* name, int size);
  void init_chars (const char* chars, const char* comp = "");  // nondegenerate chars
  void init_degen_chars (const char* chars);  // degenerate chars
  void init_degen_complement();
  void init_display (const sstring& chars);
  void set_degen_char (char degen_char, const Symbol_weight_map& p);
  void set_degen_char (char degen_char, const char* possibilities);
  void set_wildcard (char wildcard_char);
  void set_unknown (char unknown_char);

  // method to initialize a "hidden" alphabet from a "visible" alphabet
  // suppose the visible alphabet has A [nondegenerate] symbols, and the hidden_classes variable is equal to C.
  // then the hidden alphabet will have C*A symbols,
  // with each symbol of the visible alphabet representing a C-fold degenerate symbol in the hidden alphabet.
  // Symbol N in the visible alphabet is a degenerate char for symbols {N,A+N,2*A+N...(C-1)*A+N} in the hidden alphabet.
  void init_hidden (const Alphabet& alph, int hidden_classes = 1);
  void init_hidden (const Alphabet& base_alphabet, const vector<sstring>& hidden_class_labels);

  // constructor
  Alphabet (const char* name, int size);

  // methods
  int size() const { return alphabet_size; }

  // chars() & chars_toupper() methods return C strings containing nondegenerate symbols (e.g. ACGT for DNA)
  const char* nondegenerate_chars() const { return chars_lc.c_str(); }
  const char* nondegenerate_chars_toupper() const { return chars_uc.c_str(); }
  const char* degenerate_chars() const { return degen_lc.c_str(); }

  // tokens() method is an accessor for tok vector
  const vector<sstring>& tokens() const { return tok; }

  // basic query/conversion methods

  inline bool contains_strict (char c) const  // forbids degenerate chars
  {
    const int i = char2int_table[c];
    return i >= 0 && i < alphabet_size;
  }

  inline bool contains_deg (char c) const  // allows degenerate chars
  {
    const int i = char2int_table[c];
    return i > - (int) degen_lc.size() && i < alphabet_size;
  }

  inline int char2int_strict (char c) const  // randomises degenerate/unknown chars
  {
    int i = char2int_table[c];
    if (i >= alphabet_size)
      {
	if (warn_unknown) CLOGERR << "Randomising unknown characters\n";
	((Alphabet*)this)->warn_unknown = 0;  // cast away const
	i = alphabet_size > 1 ? Rnd::rnd_int (alphabet_size) : 0;  // don't call Rnd unless necessary
      }
    else if (i < 0)  // degenerate char
      {
	const Symbol_weight_map& deg_map = int2weight(i);
	const vector<int> degen_opt = map_keys (deg_map);
	const vector<Prob> degen_prob = map_values (deg_map);
	i = degen_opt.size() > 1 ? degen_opt[Rnd::choose(degen_prob)] : degen_opt[0];  // don't call Rnd unless necessary
	if (degen_opt.size() > 1)
	  {
	    if (warn_randomise) CLOGERR << "Randomising degenerate characters\n";
	    ((Alphabet*)this)->warn_randomise = 0;  // cast away const
	  }
      }
    return i;
  }
  
  inline int char2int_deg (char c) const  // converts unknown chars to (unknown_int)
  {
    int i = char2int_table[c];
    if (i >= alphabet_size)
      {
	i = unknown_int;
	if (i >= alphabet_size)
	  {
	    i = alphabet_size > 1 ? Rnd::rnd_int (alphabet_size) : 0;  // don't call Rnd unless necessary
	    CLOGERR << "I don't think '" << c << "' matches alphabet " << name << "; randomizing to '" << int2char(i) << "'\n";
	  }
	else
	  CLOGERR << "I don't think '" << c << "' matches alphabet " << name << "; setting to '" << int2char(i) << "'\n";
      }
    return i;
  }

  inline char int2char (int i) const
  {
    if (i >= 0)
      {
	if (i >= alphabet_size) return Alphabet_nonsense_char;
	if (i >= (int) chars_lc.size())
	  return i >= (int) display_lc.size() ? Alphabet_nonsense_char : display_lc[i];
	return chars_lc[i];
      }
    if (-i-1 >= (int) degen_lc.size()) return Alphabet_nonsense_char;
    return degen_lc[-i-1];
  }

  inline char int2char_uc (int i) const
  {
    int c = int2char (i);
    if (c >= 'a' && c <= 'z') c += 'A' - 'a';
    return c;
  }

  inline char unknown_char() const { return int2char (unknown_int); }

  // fancier conversion methods

  inline const Symbol_weight_map& int2weight (int i) const { return i >= 0 ? nondegen_weight[i] : degen_weight[-i-1]; }
  inline const Symbol_weight_map& char2weight (char c) const { return int2weight (char2int_deg (c)); }
  
  int weight2int (const Symbol_weight_map& w) const;
  char weight2char (const Symbol_weight_map& w) const;

  inline const Symbol_score_map& int2score (int i) const { return i >= 0 ? nondegen_score[i] : degen_score[-i-1]; }
  inline const Symbol_score_map& char2score (char c) const { return int2score (char2int_deg (c)); }

  int score2int (const Symbol_score_map& s) const;
  char score2char (const Symbol_score_map& s) const;

  // Whole-sequence conversions.
  // Each of these methods returns the second argument as a reference
  // (not really DART's most up-to-date C++ coding practise, if such a standard exists)

  Digitized_biosequence&  seq2dsq  (const Biosequence& seq,           Digitized_biosequence& dsq) const;
  Biosequence&            dsq2seq  (const Digitized_biosequence& dsq, Biosequence& seq)           const;

  Weight_profile&         seq2weight  (const Biosequence& seq,     Weight_profile& prof) const;
  Biosequence&            weight2seq  (const Weight_profile& prof, Biosequence& seq)     const;

  Score_profile&          seq2score  (const Biosequence& seq,    Score_profile& prof) const;
  Biosequence&            score2seq  (const Score_profile& prof, Biosequence& seq) const;

  Weight_profile&         dsq2weight  (const Digitized_biosequence& dsq, Weight_profile& prof) const;
  Digitized_biosequence&  weight2dsq  (const Weight_profile& prof,       Digitized_biosequence& dsq) const;   // no implementation yet

  Score_profile&          dsq2score  (const Digitized_biosequence& dsq,  Score_profile& prof) const;
  Digitized_biosequence&  score2dsq  (const Score_profile& prof,         Digitized_biosequence& dsq) const;   // no implementation yet

  // sequence conversion methods that dynamically allocate & return a new sequence

  Digitized_biosequence  new_seq2dsq  (const Biosequence& seq) const;
  Biosequence            new_dsq2seq  (const Digitized_biosequence& dsq) const;

  Weight_profile         new_seq2weight  (const Biosequence& seq) const;
  Biosequence            new_weight2seq  (const Weight_profile& prof)     const;

  Score_profile          new_seq2score  (const Biosequence& seq) const;
  Biosequence            new_score2seq  (const Score_profile& prof) const;

  Weight_profile         new_dsq2weight  (const Digitized_biosequence& dsq) const;
  Digitized_biosequence  new_weight2dsq  (const Weight_profile& prof) const;   // no implementation yet

  Score_profile          new_dsq2score  (const Digitized_biosequence& dsq) const;
  Digitized_biosequence  new_score2dsq  (const Score_profile& prof) const;  // no implementation yet

  // wildcard stuff
  bool all_wild (const Score_profile& prof) const;     // returns TRUE if Score_profile contains nothing but wildcards
  bool all_wild (const Weight_profile& prof) const;    // returns TRUE if Weight_profile contains nothing but wildcards

  Symbol_weight_map       flat_weight_map (double weight) const;      // "*" (weight = 1) or "N" (DNA, weight = .25) or "X" (protein, weight = .05)
  Symbol_score_map        flat_score_map (int score) const;           // like flat_weight_map() but in score space

  // complement stuff
  bool                    has_complement() const { return nondegen_comp.size() != 0; }
  int                     complement (int sym) const;

  Symbol_weight_map       complement_weight (const Symbol_weight_map& m) const;
  Symbol_score_map        complement_score (const Symbol_score_map& m) const;
  char                    complement_char (char c) const;

  Weight_profile          revcomp_weight_profile (const Weight_profile& p) const;    // NB does a _reverse_ complement
  Score_profile           revcomp_score_profile (const Score_profile& p) const;
  Digitized_biosequence   revcomp_digitized_sequence (const Digitized_biosequence& dsq) const;
  Biosequence             revcomp_sequence (const Biosequence& s) const;             // preserves case

  // nmer stuff
  int nmer_token (int word_len, const Digitized_biosequence::const_iterator iter);
  int nmer_token (const char* nmer);

  // regular expression helper
  sstring degenerate_seq_pattern (const Biosequence& seq) const;  // converts IUPAC degenerate sequence to regexp, e.g. "n" -> "[acgt]" for DNA

  // method to return the affinity of an alphabet for a particular character
  // used to figure out which alphabet corresponds to an anonymous input sequence (typically "is it DNA or Protein?")
  inline double char_affinity (char c) const { return char_aff_table[c]; }
};

struct Empty_alphabet : public Alphabet
{
  Empty_alphabet() : Alphabet ("", 0) { }
};

struct Simple_alphabet : public Alphabet
{
  Simple_alphabet (const char* name, const char* chars);
};

struct DNA_alphabet_type : public Alphabet
{
  DNA_alphabet_type();
};

struct RNA_alphabet_type : public Alphabet
{
  RNA_alphabet_type();
};

struct Protein_alphabet_type : public Alphabet
{
  Protein_alphabet_type();
};

struct Roman_alphabet_type : public Alphabet
{
  Roman_alphabet_type();
};

struct Text_alphabet_type : public Alphabet
{
  Text_alphabet_type();
};

extern DNA_alphabet_type DNA_alphabet;
extern RNA_alphabet_type RNA_alphabet;
extern Protein_alphabet_type Protein_alphabet;
extern Text_alphabet_type Text_alphabet;
extern Roman_alphabet_type Roman_alphabet;
extern Simple_alphabet Dummy_alphabet;  // only has one symbol

// A Metascore is an array of scores associated with each residue in a sequence,
// representing splice site predictions, masking, or any other secondary scoring source.

typedef vector<Score> Metascore;
typedef vector<Prob>  Metaprob;   // probability-space version of a Metascore

// A Meta_profile is a base class for Named_profile, holding Score_profile & Metascore vector.
// (Metascores can be used to store scores associated with sequences, like splice site predictions or masks)
struct Meta_profile
{
  Score_profile         prof_sc;  // sequence in score profile form
  vector<Metascore>     meta_sc;  // metascores
};

// Profile_flags_enum, for specifying type conversions and the like
struct Profile_flags_enum
{
  enum Profile_flags { NONE = 0, SEQ = 1, SCORE = 2, WEIGHT = 4, DSQ = 8, ALL = 15 };
};

// (name, start, end) tuple
struct NSE
{
  // data
  sstring seqname;
  int     start;
  int     end;
  // constructors
  NSE();
  NSE (const sstring& name, int start, int end);
};

// A Named_profile holds name, "cruft" (i.e. the text following the name in FASTA files),
// sequence data (in text, digitized, score and/or probability form) and Metascores.
// 
//
// No attempt is made to maintain consistency between alternative forms of sequence -- that is the caller's job.
// Sometimes, when non-IUPAC degenerate score profiles are used, it may not even be possible to interconvert.
//
// ...yes I *know* this is messy. it keeps me awake at nights. maybe one day i'll fix it. (chyeah)
//
struct Named_profile : Meta_profile, Profile_flags_enum
{
  sstring               name;     // name
  sstring               cruft;    // "cruft" (text following name in FASTA definition line)
  Biosequence           seq;      // sequence in text form
  Weight_profile        prof_w;   // sequence in probability profile form
  Digitized_biosequence dsq;      // sequence in digitised form

  // regexps for FASTA-format input
  static Regexp name_regexp;
  static Regexp quick_name_regexp;

  // regexp for Stockholm "NAME/START-END"
  static Regexp nse_regexp;

  // methods
  // clear method
  void clear();

  // FASTA IO
  void write_FASTA (ostream& out, int column_width = 50) const;
  void read_FASTA (istream& in);

  // metascore IO
  void write_metaprob (int meta_idx, ostream& out, int column_width = 50) const;
  void read_metaprob (int meta_idx, istream& in, bool expect_name = 1);

  // method to check a file for FASTA form (initial ">")
  static bool detect_FASTA (istream& in);

  // quick FASTA output wrapper
  void dump() const { write_FASTA (cout, 50); }

  // all-purpose conversion method
  void seq_update (const Alphabet& a, Profile_flags convert = ALL);

  // conversion methods
  void seq2score (const Alphabet& a) { a.seq2score (seq, prof_sc); }
  void seq2weight (const Alphabet& a) { a.seq2weight (seq, prof_w); }
  void seq2dsq (const Alphabet& a) { a.seq2dsq (seq, dsq); }

  void score2seq (const Alphabet& a) { a.score2seq (prof_sc, seq); }
  void score2weight() { prof_w = prof_sc.score2weight(); }
  void score2dsq (const Alphabet& a);    // no implementation yet

  void weight2seq (const Alphabet& a) { a.weight2seq (prof_w, seq); }
  void weight2score() { prof_sc = Score_profile (prof_w); }
  void weight2dsq (const Alphabet& a);   // no implementation yet

  void dsq2seq (const Alphabet& a) { a.dsq2seq (dsq, seq); }
  void dsq2score (const Alphabet& a) { a.dsq2score (dsq, prof_sc); }
  void dsq2weight (const Alphabet& a);   // no implementation yet

  // consensus conversion methods
  void consensus_score2seq (const Alphabet& a) { a.dsq2seq (prof_sc.consensus_dsq(), seq); }
  void consensus_weight2seq (const Alphabet& a) { a.dsq2seq (prof_w.consensus_dsq(), seq); }

  void consensus_score2dsq (const Alphabet& a) { dsq = prof_sc.consensus_dsq(); }
  void consensus_weight2dsq (const Alphabet& a) { dsq = prof_w.consensus_dsq(); }

  // reverse complement method
  Named_profile revcomp (const Alphabet& a, const char* name_suffix = "", const char* cruft_suffix = "") const;

  // metascore reset & update methods
  void reset_metascores (int n_metascores, Score default_sc = 0);
  void add_metascores (int metascore_idx, const Metascore& delta);

  // size method
  int size() const;  // does a simple consistency check between profile lengths

  // subseq extraction method
  Named_profile subseq (int start, int len) const;

  // Stockholm "NAME/START-END" parsing method
  bool name_is_NSE() const;
  NSE parse_name() const;

  // the following method, assert_profiles_consistent(), checks the various alternate profile sizes for self-consistency.
  // this is not a perfect check but it's close enough for National Lab work.
  // each profile is checked if its length is nonzero (as with size()) OR if its test flag is set.
  void assert_profiles_consistent (Profile_flags test = ALL) const;
};


// A Sequence_database is a very basic container for Named_profile's.
// NB it is generally not a great idea to instantiate a Sequence_database anywhere except
// at the top level of your program (i.e. don't give objects their own Sequence_database's)
// due to the ill-thought-out way that Alignments point to Score_profiles
// (the Alignment's pointers can be invalidated if the Sequence_database is copied).
// A general DART principle is that Sequence_database's are considered persistent,
// and iterators pointing to the Named_profile's inside will remain valid.
// This is why Sequence_database inherits from list<>, rather than e.g. vector<>
// (which tends to invalidate its iterators whenever a resizing is necessary).
//
struct Sequence_database : list<Named_profile>, Profile_flags_enum
{
  // virtual destructor
  virtual ~Sequence_database() { }

  // methods
  void write_FASTA (ostream& out, int column_width = 50) const;
  void read_FASTA (istream& in);

  void write_metaprob (int meta_idx, ostream& out, int column_width = 50) const;
  void read_metaprob (int meta_idx, istream& in);

  void to_lower();
  void to_upper();

  void seqs2scores (const Alphabet& a);
  void seqs2weights (const Alphabet& a);
  void seqs2dsqs (const Alphabet& a);

  void seqs_update (const Alphabet& a, Profile_flags convert = ALL);

  void scores2seqs (const Alphabet& a);
  void scores2weights();
  void scores2dsqs();      // no implementation yet

  void weights2seqs (const Alphabet& a);
  void weights2scores();
  void weights2dsqs();     // no implementation yet

  void dsqs2seqs (const Alphabet& a);   // no implementation yet
  void dsqs2scores();                   // no implementation yet
  void dsqs2weights();                  // no implementation yet

  void consensus_scores2seqs (const Alphabet& a);
  void consensus_weights2seqs (const Alphabet& a);

  void consensus_scores2dsqs (const Alphabet& a);       // no implementation yet
  void consensus_weights2dsqs (const Alphabet& a);      // no implementation yet

  void clear_seqs() { for_contents (list<Named_profile>, *this, prof) (*prof).seq.clear(); }
  void clear_scores() { for_contents (list<Named_profile>, *this, prof) (*prof).prof_sc.clear(); }
  void clear_weights() { for_contents (list<Named_profile>, *this, prof) (*prof).prof_w.clear(); }
  void clear_dsqs() { for_contents (list<Named_profile>, *this, prof) (*prof).dsq.clear(); }

  void reset_metascores (int n_metascores, Score default_sc = 0)
  { for_contents (list<Named_profile>, *this, prof) (*prof).reset_metascores (n_metascores, default_sc); }

  Sequence_database revcomp (const Alphabet& a, const char* name_suffix = "", const char* cruft_suffix = "") const;

  bool            matches_alphabet (const Alphabet& alphabet, double threshold = .95) const;   // returns TRUE if fraction of symbols represented by alphabet is greater than threshold
  double alphabet_similarity (const Alphabet& alphabet) const;  // each nondegenerate symbol score 1, each degenerate scores 0.5
  const Alphabet& detect_alphabet() const;  // returns DNA_alphabet, RNA_alphabet, Protein_alphabet or Text_alphabet, or dies

  void count_symbol_frequencies (vector<double>& composition) const;   // requires Digitized_biosequence's to be pre-computed by e.g. calling seqs2dsqs()
  vector<Prob> get_null_model (int alphabet_size, double pseudocount = 0.) const;
  Prob get_null_extend (double pseudocount = 0.) const;  // optimal extension prob. for null model
  double mean_length() const;
  int total_residues() const;

  // helper to append an entire Sequence_database
  // returns iterator pointing to first new sequence
  virtual iterator append_sequence_database (const Sequence_database& seq_db);
};

// a Sequence_database_index is a by-name indexing of a Sequence_database
//
struct Sequence_database_index
{
  vector<Named_profile*> profile;
  vector<sstring>        name;
  Phonebook              name2profile_index;

  Sequence_database_index() { }
  Sequence_database_index (Sequence_database& seq_db) { update (seq_db); }

  void update (Sequence_database& seq_db);  // rebuilds profile, name & name2profile_index
  int  size() const { return profile.size(); }
  void assert_names_unique (bool test_NSE_conflict = TRUE) const;

  Named_profile* name2profile (const sstring& name) const;
};

// FASTA_sequence_database is a convenient class for reading in FASTA files

class FASTA_sequence_database : public Sequence_database
{
  // data
private:
  const Alphabet* my_alphabet;
public:
  Sequence_database_index index;

  // accessors
  const Alphabet& alphabet() const { return *my_alphabet; }
  int size() const { return index.size(); }    // override size method for speed (vector::size() is constant)

  const Named_profile& get_seq (int i) const { return *index.profile[i]; }

  // constructor
  FASTA_sequence_database (const char* filename = 0,
			   const Alphabet* alphabet = 0,
			   Profile_flags convert = ALL);   // calls read_from_file()

  // I/O
  void read_FASTA (istream& in);  // calls superclass method & then updates index, but does not call seqs_update()
  void read_from_file (const char* filename,
		       const Alphabet* alphabet = 0,
		       Profile_flags convert = ALL);   // calls seqs_update()

  // helpers for updating index, alphabet
  void update_index() { index.update (*this); }
  void update_alphabet (const Alphabet* alphabet = 0);

  // helper for manually adding sequences
  void add_seq (const char* name, const char* seq, Profile_flags convert = ALL);  // calls seq_update() on new sequence

  // helper to append an entire Sequence_database (overrides Sequence_database method)
  // returns iterator pointing to first new sequence
  iterator append_sequence_database (const Sequence_database& seq_db);
};

#endif
