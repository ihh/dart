#include "ecfg/ecfgsexpr.h"
#include "util/unixenv.h"
#include "hmm/transmat.h"
#include <iostream>
#include <ext/hash_map>
#include <cstring>
#include <vector>
// path to default chain file, from dart root
#define DEFAULT_CHAIN_FILE "data/handalign/prot1.hsm"

// wrapper for Alphabet::char2int
int char2int (const Alphabet& alphabet, const sstring& s, const char* desc)
{
  if (s.size() < 1)
    THROWEXPR (desc << " state string '" << s << "' has too few characters (must be exactly one char)");

  if (s.size() > 1)
    THROWEXPR (desc << " state string '" << s << "' has too many characters (must be exactly one char)");

  const char c = s[0];
  if (!alphabet.contains_strict (c))
    THROWEXPR (desc << " character '" << c << "' not found in alphabet");

  return alphabet.char2int_strict (c);
}

using namespace std;


struct eqstr{
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1,s2)==0;
  }
};

struct eqvec{
  bool operator()(const vector<int> v1, const vector<int> v2) const {
    return (v1[0] == v2[0] && 
	    v1[1] == v2[1] && 
	    v1[2] == v2[2] && 
	    v1[3] == v2[3] && 
	    v1[4] == v2[4] );
      }
};


class example: public vector<int>
{
public:
  uint8_t name_[8];
};

namespace __gnu_cxx {
  template<> struct hash<example > {
    size_t operator()(const example& k) const {
      size_t hashval = 0;
      for (int i = 0; i < sizeof(k.name_); ++i) {
	hashval = 5 * hashval + k.name_[i];
      }
      return hashval;
    }
  };
}  // namespace __gnu_cxx


// main program
int main (int argc, char** argv)
{
//   // initialise the options parser
//   INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
// 		  "read a substitution matrix and exponentiate it\n");

//   sstring default_chain_filename;
//   default_chain_filename << Dart_Unix::get_DARTDIR() << '/' << DEFAULT_CHAIN_FILE;

//   sstring chain_filename, src_char, dest_char;
//   double time;

//   opts.add ("c -chain-filename", chain_filename = default_chain_filename, "chain file");
//   opts.add ("s -source-state", src_char = "a", "source state character");
//   opts.add ("d -dest-state", dest_char = "a", "destination state character");
//   opts.add ("t -time", time = 1., "rate matrix multiplier");

//   // parse the command line & do stuff
//   try
//     {
//       // parse command line
//       if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
//     }
//   catch (const Dart_exception& e)
//     {
//       cerr << opts.short_help();
//       cerr << e.what();
//       exit(1);
//     }

//   // do stuff
//   try
//     {
//       // read chain file, get SExpr
//       SExpr_file ecfg_sexpr_file (chain_filename.c_str());
//       SExpr& ecfg_sexpr = ecfg_sexpr_file.sexpr;

//       // init alphabet, rate matrix
//       Alphabet alphabet ("uninitialized", 1);  // dummy init of alphabet (from seq/biosequence.h)
//       Irrev_EM_matrix rate_matrix (1, 1);  // dummy init of rate_matrix (from irrev/irrev_em_matrix.h)
//       ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, ecfg_sexpr);

// //       // print alphabet tokens, in order
//       cout << "Alphabet symbols: " << alphabet.tokens() << '\n';


// //       // get equilibrium distribution
//        const vector<double> eqm_prob = rate_matrix.create_prior();

//       // print equilibrium distribution
//       cout << "Equilibrium distribution: " << eqm_prob << '\n';

//       // convert alphabet symbols to state indices
// 	  const int src_state = 0;
// 	  const int dest_state = 1;

// 	  //       // exponentiate matrix
// 	  const array2d<double> cond_submat = rate_matrix.create_conditional_substitution_matrix (time);
	  
	  
	  
// 	  //       // extract required element
// 	  const double result = cond_submat (src_state, dest_state);
	  
// 	  cout << result<<endl;
//       // print required element
// 	  //      cout << "exp(R*" << time << ")_{" << src_char << ',' << dest_char << "} = " << result << '\n';
//     }
//   catch (const Dart_exception& e)
//     {
//       cerr << e.what();
//       exit(1);
//     }

//   array2d<double> test_matrix(3,3); 
//   test_matrix.fill(.3333); 
//   std::cout<<"The transition matrix: \n"; 
//   test_matrix.write_rowcol(std::cout);

//   Transition_probs t(3); 
//   t.transition(0,1) = .5; 
//   t.transition(1,2) = .5; 

//   //  t.assign_transition_matrix(test_matrix); 
//   t.print_dotfile_edges(std::cout); 
  
//   vector<int> nulls; 
//   nulls.push_back(1); 
//   Transition_methods m; 
//   Concrete_transition_probs ct = m.eliminate(t, nulls); 

//   ct.print_dotfile_edges(std::cout); 
//   std::cout<<ct.transition(0,1) << endl; 
  using namespace __gnu_cxx; 
  typedef __gnu_cxx::hash_map<example,double> Example1HashType;
  hash_map<const char*, int, hash<const char*>, eqstr> months;
  months["january"] = 31;
  months["february"] = 28;
  months["march"] = 31;
  months["april"] = 30;
  months["may"] = 31;
  months["june"] = 30;
  months["july"] = 31;
  months["august"] = 31;
  months["september"] = 30;
  months["october"] = 31;
  months["november"] = 30;
  months["december"] = 31;
  cout << "september -> " << months["september"] << endl;
  cout << "april     -> " << months["april"] << endl;
  cout << "june      -> " << months["june"] << endl;
  cout << "november  -> " << months["november"] << endl; 
  cout << "march     -> " << months["march"] << endl; 
  
  Example1HashType DP;
  
  example m;
  m.push_back(3);     m.push_back(3);     m.push_back(3);     m.push_back(3);     m.push_back(3); 
  DP[m] = 1.01;
  cout<< "value: " << DP[m] <<endl; 
  
  return 0;
}
