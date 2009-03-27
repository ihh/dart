#include <ostream>
#include "handel/ehmm.h"
#include "util/vector_output.h"

// Default TAlphabet
struct DefaultTAlphabet : TAlphabet
{
  DefaultTAlphabet (int terminals) : TAlphabet (terminals, '?')
  {
    for (int x = 0; x < terminals; ++x)
      (*this)[x] = 'A' + x;
  }
};

// Function to test creation of an ETree from a string
bool testTreeString (const char* treestr, bool leafnames = TRUE, int nodes = -1, const char* expect = 0)
{
  bool ok = TRUE;
  const basic_string<char> initString (treestr);
  const ETree tree = ETree::string2etree (initString);
  cerr << " [ETree initialised from \"" << initString << "\"";

  // test validity
  if (!tree.valid())
    {
      cerr << " is invalid";
      ok = FALSE;
    }
  else
    cerr << " is valid";
  
  // test number of nodes
  const int treenodes = tree.nodes();
  cerr << ", has " << treenodes << " nodes";
  if (nodes >= 0)
    if (treenodes != nodes)
      {
	cerr << " where " << nodes << " expected";
	ok = FALSE;
      }
    else
      cerr << " as expected";
  else
    cerr << " assumed ok";

  // test string
  const basic_string<char> treeString = tree.etree2string (leafnames);
  cerr << ", string \"" << treeString << "\"";
  if (expect)
    {
      const basic_string<char> expectString (expect);
      if (treeString != expectString)
	{
	  cerr << " where \"" << expect << "\" expected";
	  ok = FALSE;
	}
      else
	cerr << " expected";
    }
  else if (treeString != initString)
    {
      cerr << " incorrect";
      ok = FALSE;
    }
  else
    cerr << " correct";

  // return
  if (ok)
    cerr << ", ok]\n";
  else
    cerr << ", NOT OK]\n";
  return ok;
}

// ETree subclass, created from a string
struct ETreeFromString : ETree
{
  ETreeFromString (const char* init) : ETree (string2etree (basic_string<char> (init))) { }
};

// Function to print an EStateList
basic_string<char> estatelist2string (const EStateList& estateList, const TAlphabet& alph)
{
  basic_string<char> s;
  for_const_contents (EStateList, estateList, estate)
    {
      if (s.size()) s.push_back (' ');
      s.append (estate->estate2string (alph));
    }
  return s;
}

// Function to print an EStateIterator & count the elements
basic_string<char> estateiter2string (EStateIterator& estateIter, int nodes, const TAlphabet& alph, int* count = 0)
{
  if (count)
    *count = 0;
  basic_string<char> s;
  EState estate (nodes);
  for (estateIter.reset(); !estateIter.finished(); estateIter.next())
    {
      estateIter.getCurrent (estate);
      if (s.size()) s.push_back (' ');
      s.append (estate.estate2string (alph));
      if (count)
	++*count;
    }
  return s;
}

// Function to dump a CollapsedESpaceSubsets object to standard error
void dumpCollapsedESubsets (const TSpace& tspace, const ETree& tree, const ESubtrees& subtrees, CollapsedESpaceSubsets& subsets)
{
  // make DefaultTAlphabet
  const DefaultTAlphabet T (tspace.terminals);

  // make refs of A, B
  vector<CollapsedESpaceSubsets::ABMap>& A = subsets.A;
  vector<CollapsedESpaceSubsets::ABMap>& B = subsets.B;
      
  // show A and B sets
  for (ENode n = 0; n < tree.nodes(); ++n)
    {
      // show node & parent
      cerr << "   Node " << n << ", parent " << tree.parent[n];
      cerr << ", children {" << subtrees.children[n] << "}";
      cerr << ", descendants {" << subtrees.descendants[n] << "}\n";

      // show A
      // Start, End, Wait
      cerr << "    A[" << n << "][Start]: {" << estatelist2string (A[n][TSpaceEnum::TransStart], T) << "}\n";
      cerr << "    A[" << n << "][End]: {" << estatelist2string (A[n][TSpaceEnum::TransEnd], T) << "}\n";
      cerr << "    A[" << n << "][Wait]: {" << estatelist2string (A[n][TSpaceEnum::TransWait], T) << "}\n";

      // Match
      for (TTerm x = 0; x < tspace.terminals; ++x)
	for (TTerm y = 0; y < tspace.terminals; ++y)
	  {
	    cerr << "    A[" << n << "][Match";
	    if (tspace.terminals > 1)
	      cerr << "(" << x << "," << y << ")";
	    cerr << "]: {" << estatelist2string (A[n][tspace.TransMatchXY(x,y)], T) << "}\n";
	  }

      // Delete
      for (TTerm x = 0; x < tspace.terminals; ++x)
	{
	  cerr << "    A[" << n << "][Delete";
	  if (tspace.terminals > 1)
	    cerr << "(" << x << ")";
	  cerr << "]: {" << estatelist2string (A[n][tspace.TransDeleteX(x)], T) << "}\n";
	}


      // show B
      // Start
      cerr << "    B[" << n << "][Start]: {" << estatelist2string (B[n][TSpaceEnum::TransStart], T) << "}\n";

      // Match
      for (TTerm x = 0; x < tspace.terminals; ++x)
	for (TTerm y = 0; y < tspace.terminals; ++y)
	  {
	    cerr << "    B[" << n << "][Match";
	    if (tspace.terminals > 0)
	      cerr << "(" << x << "," << y << ")";
	    cerr << "]: {" << estatelist2string (B[n][tspace.TransMatchXY(x,y)], T) << "}\n";
	  }

      // Insert
      for (TTerm y = 0; y < tspace.terminals; ++y)
	{
	  cerr << "    B[" << n << "][Insert";
	  if (tspace.terminals > 0)
	    cerr << "(" << y << ")";
	  cerr << "]: {" << estatelist2string (B[n][tspace.TransInsertY(y)], T) << "}\n";
	}
    }
}

// Function to dump a CollapsedEStateIterators object to standard error
void dumpCollapsedEStateIterators (const TSpace& tspace, const ETree& tree, const ESubtrees& subtrees, CollapsedEStateIterators& iters)
{
  // make DefaultTAlphabet
  const DefaultTAlphabet T (tspace.terminals);

  // make shorthand for tree nodes
  const int N = tree.nodes();

  // show A and B sets
  for (ENode n = 0; n < tree.nodes(); ++n)
    {
      // show node & parent
      cerr << "   Node " << n << ", parent " << tree.parent[n];
      cerr << ", children {" << subtrees.children[n] << "}";
      cerr << ", descendants {" << subtrees.descendants[n] << "}\n";

      // show A
      // Start, End, Wait
      cerr << "    A[" << n << "][Start]: {" << estateiter2string (iters.enodeIter[n].tstateIter[TSpaceEnum::TransStart].A, N, T) << "}\n";
      cerr << "    A[" << n << "][End]: {" << estateiter2string (iters.enodeIter[n].tstateIter[TSpaceEnum::TransEnd].A, N, T) << "}\n";
      cerr << "    A[" << n << "][Wait]: {" << estateiter2string (iters.enodeIter[n].tstateIter[TSpaceEnum::TransWait].A, N, T) << "}\n";

      // Match
      for (TTerm x = 0; x < tspace.terminals; ++x)
	for (TTerm y = 0; y < tspace.terminals; ++y)
	  {
	    cerr << "    A[" << n << "][Match";
	    if (tspace.terminals > 1)
	      cerr << "(" << x << "," << y << ")";
	    cerr << "]: {" << estateiter2string (iters.enodeIter[n].tstateIter[tspace.TransMatchXY(x,y)].A, N, T) << "}\n";
	  }

      // Delete
      for (TTerm x = 0; x < tspace.terminals; ++x)
	{
	  cerr << "    A[" << n << "][Delete";
	  if (tspace.terminals > 1)
	    cerr << "(" << x << ")";
	  cerr << "]: {" << estateiter2string (iters.enodeIter[n].tstateIter[tspace.TransDeleteX(x)].A, N, T) << "}\n";
	}


      // show B
      // Start
      cerr << "    B[" << n << "][Start]: {" << estateiter2string (iters.enodeIter[n].tstateIter[TSpaceEnum::TransStart].B, N, T) << "}\n";

      // Match
      for (TTerm x = 0; x < tspace.terminals; ++x)
	for (TTerm y = 0; y < tspace.terminals; ++y)
	  {
	    cerr << "    B[" << n << "][Match";
	    if (tspace.terminals > 0)
	      cerr << "(" << x << "," << y << ")";
	    cerr << "]: {" << estateiter2string (iters.enodeIter[n].tstateIter[tspace.TransMatchXY(x,y)].B, N, T) << "}\n";
	  }

      // Insert
      for (TTerm y = 0; y < tspace.terminals; ++y)
	{
	  cerr << "    B[" << n << "][Insert";
	  if (tspace.terminals > 0)
	    cerr << "(" << y << ")";
	  cerr << "]: {" << estateiter2string (iters.enodeIter[n].tstateIter[tspace.TransInsertY(y)].B, N, T) << "}\n";
	}

      // show Amd, Bmi, AmdBmi
      for (TTerm x = 0; x < tspace.terminals; ++x)
	{
	  cerr << "    Amd[" << n << "]";
	  if (tspace.terminals > 0)
	    cerr << "[" << x << "]";
	  cerr << ": {" << estateiter2string (iters.enodeIter[n].ttermIter[x].Amd, N, T) << "}\n";
	  cerr << "    Bmi[" << n << "]";
	  if (tspace.terminals > 0)
	    cerr << "[" << x << "]";
	  cerr << ": {" << estateiter2string (iters.enodeIter[n].ttermIter[x].Bmi, N, T) << "}\n";
	  cerr << "    AmdBmi[" << n << "]";
	  if (tspace.terminals > 0)
	    cerr << "[" << x << "]";
	  cerr << ": {" << estateiter2string (iters.enodeIter[n].ttermIter[x].AmdBmi, N, T) << "}\n";
	}

      // show Bsi
      cerr << "    Bsi[" << n << "]: {" << estateiter2string (iters.enodeIter[n].Bsi, N, T) << "}\n";
    }
}

// Function to build collapsed EHMM state spaces
bool makeCollapsedEHMM (const TSpace& tspace, const ETree& tree, int expectedConditionalStates, int expectedJointStates,
			bool showCondmatrix = FALSE, bool showJointmatrix = FALSE)
{
  // set ok flag
  bool ok = TRUE;

  // make DefaultTAlphabet
  const DefaultTAlphabet T (tspace.terminals);

  // make ESubtrees
  const ESubtrees subtrees (tree);

  // make CollapsedESpaceSubsets
  cerr << " [Building EHMM for tree " << tree.etree2string();
  cerr << " with " << tree.nodes() << " nodes, leaves {" << subtrees.leaves << "}, ";
  cerr << tspace.terminals << " terminals]\n";
  CollapsedESpaceSubsets subsets (tspace, tree);

  // show joint state space
  const EStateList joint = subsets.joint (tspace);
  const int js = joint.size();
  cerr << "  [" << js << " joint states";
  if (js != expectedJointStates) { cerr << " (expected " << expectedJointStates << ")"; ok = FALSE; }
  cerr << ": " << estatelist2string (joint, T) << "]\n";

  // show conditional state space
  const EStateList conditional = subsets.conditional (tspace);
  const int cs = conditional.size();
  cerr << "  [" << cs << " conditional states";
  if (cs != expectedConditionalStates) { cerr << " (expected " << expectedConditionalStates << ")"; ok = FALSE; }
  cerr << ": " << estatelist2string (conditional, T) << "]\n";

  // if not ok, dump subsets
  if (!ok)
    dumpCollapsedESubsets (tspace, tree, subtrees, subsets);

  // make iterators
  CollapsedEStateIterators iters (tspace, tree);
  EStateUnionIterator condIter = iters.conditional();
  EStateUnionIterator jointIter = iters.joint();
  int condSize, jointSize;

  // show joint state space
  const basic_string<char> jstr = estateiter2string (jointIter, tree.nodes(), T, &jointSize);
  cerr << "  [joint iteration (" << jointSize << " states";
  if (jointSize != expectedJointStates) { cerr << ", expected " << expectedJointStates; ok = FALSE; }
  cerr << "): " << jstr << "]\n";

  // repeat the iteration, to check reset() works
  const basic_string<char> jstr2 = estateiter2string (jointIter, tree.nodes(), T);
  if (jstr2 != jstr)
    {
      ok = FALSE;
      cerr << "  [second joint iteration different from first: " << jstr2 << "]\n";
    }

  // show conditional state space
  const basic_string<char> cstr = estateiter2string (condIter, tree.nodes(), T, &condSize);
  cerr << "  [conditional iteration (" << condSize << " states";
  if (condSize != expectedConditionalStates) { cerr << ", expected " << expectedConditionalStates; ok = FALSE; }
  cerr << "): " << cstr << "]\n";

  // dump iterator tree if not ok
  if (!ok)
    dumpCollapsedEStateIterators (tspace, tree, subtrees, iters);

  // make conditional and joint matrices
  const CollapsedEMatrix condmat (conditional, tspace, tree);
  const CollapsedEMatrix jointmat (joint, tspace, tree);

  // show conditional and/or joint matrices
  if (ok && showCondmatrix)
    {
      cerr << "  Conditional matrix:\n";
      for_const_contents (vector<CollapsedETrans>, condmat.etrans, etrans)
	cerr << "   " << etrans->etrans2string (T) << "\n";
    }

  if (ok && showJointmatrix)
    {
      cerr << "  Joint matrix:\n";
      for_const_contents (vector<CollapsedETrans>, jointmat.etrans, etrans)
	cerr << "   " << etrans->etrans2string (T) << "\n";
    }

  // return
  return ok;
}

// Main
int main (int argc, char** argv)
{
  bool ok = TRUE;

  // test ETrees
  cerr << "[testing ETree initialisation]\n";
  ok = testTreeString ("(A)", TRUE, 1, "(A);") && ok;
  ok = testTreeString ("();", FALSE, 1) && ok;
  ok = testTreeString ("((A));", TRUE, 2) && ok;
  ok = testTreeString ("(());", FALSE, 2) && ok;
  ok = testTreeString ("((A,B));", TRUE, 3) && ok;
  ok = testTreeString ("((,));", FALSE, 3) && ok;
  ok = testTreeString ("((((antU))))", FALSE, 4, "(((())));") && ok;
  ok = testTreeString ("(((A,B,C),D));", TRUE, 6) && ok;
  ok = testTreeString ("(((,,),));", FALSE, 6) && ok;

  // create a TSpace with one terminal
  const TSpace oneTerm;

  // try CollapsedESpaceSubsets for one-branch tree (root->subroot)
  const ETreeFromString oneBranch ("()");
  ok = makeCollapsedEHMM (oneTerm, oneBranch, 6, 4, TRUE, FALSE) && ok;

  // try CollapsedESpaceSubsets for two-branch tree
  const ETreeFromString twoBranch ("(())");
  ok = makeCollapsedEHMM (oneTerm, twoBranch, 11, 7, TRUE, TRUE) && ok;

  // try CollapsedESpaceSubsets for three-branch forking tree
  const ETreeFromString threeBranchFork ("((,))");
  ok = makeCollapsedEHMM (oneTerm, threeBranchFork, 20, 12, FALSE, TRUE) && ok;

  // try CollapsedESpaceSubsets for three-branch linear tree
  const ETreeFromString threeBranchLinear ("((()))");
  ok = makeCollapsedEHMM (oneTerm, threeBranchLinear, 22, 14) && ok;

  // try CollapsedESpaceSubsets for four-branch star tree
  const ETreeFromString fourBranchStar ("((,,))");
  ok = makeCollapsedEHMM (oneTerm, fourBranchStar, 37, 21) && ok;

  // try CollapsedESpaceSubsets for four-branch forking tree
  const ETreeFromString fourBranchFork ("((,),)");
  ok = makeCollapsedEHMM (oneTerm, fourBranchFork, 20, 12, FALSE, TRUE) && ok;

  // create a TSpace with two terminals
  const TSpace twoTerm (2);

  // try CollapsedESpaceSubsets for one-branch tree (root->subroot)
  ok = makeCollapsedEHMM (twoTerm, oneBranch, 11, 5, FALSE, FALSE) && ok;

  // try CollapsedESpaceSubsets for two-branch tree
  ok = makeCollapsedEHMM (twoTerm, twoBranch, 37, 15, FALSE, FALSE) && ok;

  // print ok message
  if (ok)
    cout << "ok\n";
  else
    cout << "not ok\n";
  
  return 0;
}
