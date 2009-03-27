#include "scfg/foldenv.h"
#include "util/logfile.h"
#include "util/Regexp.h"

// global variables

Fold_envelope env, env2;

// test functions

sstring test_pseudovecs (Fold_envelope& env)
{
  sstring errstr;
  // commented out because pseudovecs have now been swapped in for vectors - ihh, 1/24/2005
  // (NB this function's a bit screwed anyway, as the vectors are reversed w.r.t. the pseudovecs)
  /*
  for (int subseq_idx = 0; subseq_idx < env.subseqs(); ++subseq_idx)
    {
      const Subseq& subseq = env.subseq[subseq_idx];
      int i;

      Subseq::Bifurc_in_pseudovec bif_in;
      bif_in.init (env, subseq_idx);
      if (bif_in.size() != subseq.bif_in.size()) {
	errstr << "bif_in: Pseudovec size = " << bif_in.size() << ", vector size = " << subseq.bif_in.size() << "\n";
	goto HELL;
      }
      i = subseq.bif_in.size() - 1;
      for_const_contents (Subseq::Bifurc_in_pseudovec, bif_in, bif_in_iter)
	{
	  Subseq::Bifurc_in b1 = *bif_in_iter, b2 = subseq.bif_in[i--];
	  if (b1.l != b2.l || b1.r != b2.r) {
	    errstr << "bif_in: subseq=" << subseq_idx  << ", b1=(" << b1.l << "," << b1.r << "), b2=(" << b2.l << "," << b2.r << ")\n";
	    goto HELL;
	  }
	  // else cerr << " [bif_in: subseq=" << subseq_idx << ", b1=(" << b1.l << "," << b1.r << "), b2=(" << b2.l << "," << b2.r << "); good]\n";
	}

      Subseq::Bifurc_outl_pseudovec bif_outl;
      bif_outl.init (env, subseq_idx);
      if (bif_outl.size() != subseq.bif_out_l.size()) {
	errstr << "bif_outl: Pseudovec size = " << bif_outl.size() << ", vector size = " << subseq.bif_out_l.size() << "\n";
	goto HELL;
      }
      i = 0;
      for_const_contents (Subseq::Bifurc_outl_pseudovec, bif_outl, bif_outl_iter)
	{
	  Subseq::Bifurc_out_l b1 = *bif_outl_iter, b2 = subseq.bif_out_l[i++];
	  if (b1.out != b2.out || b1.l != b2.l) {
	    errstr << "bif_outl: subseq=" << subseq_idx  << ", b1=(" << b1.out << "," << b1.l << "), b2=(" << b2.out << "," << b2.l << ")\n";
	    goto HELL;
	  }
	  // else cerr << " [bif_outl: subseq=" << subseq_idx << ", b1=(" << b1.out << "," << b1.l << "), b2=(" << b2.out << "," << b2.l << "); good]\n";
	}

      Subseq::Bifurc_outr_pseudovec bif_outr;
      bif_outr.init (env, subseq_idx);
      if (bif_outr.size() != subseq.bif_out_r.size()) {
	errstr << "bif_outr: Pseudovec size = " << bif_outr.size() << ", vector size = " << subseq.bif_out_r.size() << "\n";
	goto HELL;
      }
      i = 0;
      for_const_contents (Subseq::Bifurc_outr_pseudovec, bif_outr, bif_outr_iter)
	{
	  Subseq::Bifurc_out_r b1 = *bif_outr_iter, b2 = subseq.bif_out_r[i++];
	  if (b1.out != b2.out || b1.r != b2.r) {
	    errstr << "bif_outr: subseq=" << subseq_idx  << ", b1=(" << b1.out << "," << b1.r << "), b2=(" << b2.out << "," << b2.r << ")\n";
	    goto HELL;
	  }
	  // else cerr << " [bif_outr: subseq=" << subseq_idx << ", b1=(" << b1.out << "," << b1.r << "), b2=(" << b2.out << "," << b2.r << "); good]\n";
	}
    }
  cerr << " [bifurcation vectors match pseudovectors; good]\n";
 HELL:
  */
  return errstr;
}


void test_serial()
{
  const sstring env_str = env.to_envelope_string();
  env2.initialise_from_envelope_string (env_str);
  const sstring env2_str = env2.to_envelope_string();
  if (env != env2)
    {
      cerr << "Original envelope:\n";
      env.dump (cerr);
      cerr << "Serialised envelope:\n";
      env2.dump (cerr);
      cerr << "Original envelope string:\n" << env_str << "\n";
      if (env2_str != env_str)
	cerr << "Serialised envelope string:\n" << env2_str << "\n";
      else
	cerr << "(serialised envelope string is identical)\n";
      THROWEXPR ("Serialisation failed");
    }
}
void init_env (const char* s, int nsub = -1, int nbif = -1)
{
  sstring fs (s);
  // swap "<" and ">", since we changed these in the code -- ihh, 18 August 2003
  // should really fix the calls to init_env() instead -- ihh, 18 August 2003
  for_contents (sstring, fs, c)
    if (*c == '>') *c = '<'; else if (*c == '<') *c = '>';
  // do the test
  env.initialise_from_fold_string (fs);
  sstring errstr;
  if (!env.test_consistent()) errstr << "Envelope inconsistent";
  else if (!env.test_fold_string (fs)) errstr << "Envelope doesn't match fold string";
  else if (nsub >= 0 && nsub != (int) env.subseq.size()) errstr << "Expected " << nsub << " subsequences\n";
  else if (nbif >= 0 && nbif != env.bifurcations()) errstr << "Expected " << nbif << " bifurcations; found " << env.bifurcations();
  else errstr = test_pseudovecs (env);
  if (errstr.size())
    THROWEXPR ("Envelope:\n" << Fold_envelope::Dump(env) << "Fold:\n" << s << "\n" << errstr << "\n");
  cerr << "(parsed \"" << s << "\" into " << env.subseq.size() << " subseqs, " << env.bifurcations() << " bifurcations)\n";
  test_serial();
}
void init_full_env (int seqlen)
{
  // thanks to maple for these formulae (i knew them, once):
  const int nsub = (seqlen + 1) * (seqlen + 2) / 2;
  const int nbif = (seqlen * 11 + (seqlen * seqlen * seqlen)) / 6 + (seqlen * seqlen) + 1;
  env.initialise_full (seqlen);
  sstring errstr;
  if (!env.test_consistent()) errstr << "Envelope inconsistent";
  else if (nsub != (int) env.subseq.size()) errstr << "Expected " << nsub << " subsequences\n";
  else if (nbif != env.bifurcations()) errstr << "Expected " << nbif << " bifurcations; found " << env.bifurcations();
  else errstr = test_pseudovecs (env);
  if (errstr.size())
    THROWEXPR ("Full envelope for length " << seqlen << ":\n" << Fold_envelope::Dump(env) << "\n" << errstr << "\n");
  cerr << "(created length " << seqlen << " full envelope with " << env.subseq.size() << " subseqs, " << env.bifurcations() << " bifurcations)\n";
  test_serial();
}

// main

int main (int argc, char** argv)
{
  try
    {
      // first try breaking the parser

      // For some reason, this exception is not being caught. This is worrying.
      // Fuck knows why it is happening. It started after I created pairenv.[h/cc] and pairhmmem.[h/cc]
      // by moving/copying code from paircfgem.[h/cc] and foldenv.[h/cc], respectively.
      // It does not appear to be due to the nesting of one 'try...catch' clause within another.
      // A worrying prospect is that it is a memory violation somewhere else, but other tests seem to work...
      // Commenting it out for now.
      // Bugger.
      // Still, added 'testcatch.cc' in dart/src/util, which appears to reproduce the problem.
      // So that's something.
      // -- ihh, 26 September 2003

      // This was because 'gcc' was aliased to /usr/bin/gcc in ~/.cshrc, but make was not seeing this,
      // so it was using an old gcc and screwing up.
      // So, I'm putting the test back in. --ihh, 26 Sept 2003 (a bit later in the evening)

      Regexp parser_re ("Too many [<>]'s");
      
      try
	{
	  init_env (">><<<");
	  THROWEXPR ("Parser failed to break");
	}
      catch (const Dart_exception& e)
	{
	  if (!parser_re.Match (e.what())) THROWEXPR (e.what());
	  cerr << "(parser broke: that's good)\n";
	}

      try
	{
	  init_env (">...");
	  THROWEXPR ("Parser failed to break");
	}
      catch (const Dart_exception& e)
	{
	  if (!parser_re.Match (e.what())) THROWEXPR (e.what());
	  cerr << "(parser broke again: v happy)\n";
	}

      cerr << "(let's try some consistency checks)\n";
      
      init_env (">.><<");
      init_env (".>.>>>>>...<<<.><..>>..<<<<.<...");
      init_env (".>.>>>>>...<<<.><..>>..<<<<.<..><><><..>.>>>>>...<<<.><..>>..<<<<.<...");

      init_full_env (10);

      cerr << "(let's check a few envelope sizes)\n";

      init_env (">>><<<", 7 /*null*/ + 3 /*><*/, 1*7 + 2*3);
      init_env (">>><<<..>>><<<",
		15 /*null*/ + 6 /*><*/ + 2 /*.*/ + 1 /*..*/ + 2 /*><.*/ + 2 /*><..*/ + 1 /*all*/,
		1*15 + 2*6 + 2*2 + 3*1 + 3*2 + 4*2 + 5*1);

      init_env (".......", 8 + 7 + 6 + 5 + 4 + 3 + 2 + 1, 8*1 + 7*2 + 6*3 + 5*4 + 4*5 + 3*6 + 2*7 + 1*8);
      
      cerr << "(good -- that's what we expected. now what about the boundary conditions)\n";

      init_env ("");
      init_env (".");
      init_env ("><");

      init_full_env (0);
      init_full_env (1);
      init_full_env (-1);

      cerr << "(that's it -- envelopes are as slick as ever)\n";
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << "ok\n";
  return 0;
}
