#include "handel/ehmm.h"

// TSpace

TSpace::TSpace (int terminals)
  : terminals (terminals)
{ }

vector<TState> TSpace::allInsertY() const
{
  vector<TState> res;
  for (TTerm t = 0; t < terminals; ++t)
    res.push_back (TransInsertY (t));
  return res;
}

vector<TState> TSpace::allDeleteX() const
{
  vector<TState> res;
  for (TTerm t = 0; t < terminals; ++t)
    res.push_back (TransDeleteX (t));
  return res;
}

vector<TState> TSpace::allMatchXY() const
{
  vector<TState> res;
  for (TTerm xterm = 0; xterm < terminals; ++xterm)
    for (TTerm yterm = 0; yterm < terminals; ++yterm)
      res.push_back (TransMatchXY (xterm, yterm));
  return res;
}

vector<TState> TSpace::allMatchX (TTerm xterm) const
{
  vector<TState> res;
  for (TTerm yterm = 0; yterm < terminals; ++yterm)
    res.push_back (TransMatchXY (xterm, yterm));
  return res;
}

vector<TState> TSpace::allMatchY (TTerm yterm) const
{
  vector<TState> res;
  for (TTerm xterm = 0; xterm < terminals; ++xterm)
    res.push_back (TransMatchXY (xterm, yterm));
  return res;
}

// TPath

basic_string<char> TPath::tpath2string() const
{
  basic_string<char> s;
  for (int i = 0; i < (int) tstate.size(); ++i)
    {
      if (i > 0)
	s.append ("->");
      s.push_back (tstate2char (tstate[i]));
    }
  return s;
}

basic_string<char> TPath::tpath2string (const TAlphabet& alphabet) const
{
  const TSpace tspace (alphabet.size());  // create temporary TSpace
  basic_string<char> s;
  for (int i = 0; i < (int) tstate.size(); ++i)
    {
      if (i > 0)
	s.append ("->");
      s.append (tspace.tstate2string (tstate[i], alphabet));
    }
  return s;
}

// ETree

ETree::ETree (int nodes)
  : parent (nodes, 0)
{
  if (nodes > 0)
    parent[0] = -1;
}

ETree::ETree (const ETree& etree)
  : parent (etree.parent)
{ }

bool ETree::valid() const
{
  if (nodes() == 0) return FALSE;  // must have at least the subroot
  if (parent[0] != -1) return FALSE;  // subroot (uniquely) must have root as parent
  vector<ENode> lastchild (nodes(), -1);
  vector<bool> passed (nodes(), FALSE);
  for (ENode node = 1; node < nodes(); ++node)
    {
      const ENode p = parent[node];
      if (p < 0 || p >= node || passed[p])
	return FALSE;  // nodes must be sorted in preorder & depth-first, and all connected
      if (lastchild[p] >= 0) passed[lastchild[p]] = TRUE;  // pass sibling in depth-first traversal
      lastchild[p] = node;
    }
  return TRUE;
}

bool ETree::binary() const
{
  const vector<int> nchild = nchildvec();
  for_const_contents (vector<int>, nchild, n)
    if (*n != 0 && *n != 2)
      return FALSE;
  return TRUE;
}

vector<ENode> ETree::leaves() const
{
  vector<int> res;
  const vector<int> nchild = nchildvec();
  for (ENode n = 0; n < nodes(); ++n)
    if (nchild[n] == 0)
      res.push_back (n);
  return res;
}

vector<ENode> ETree::internals() const
{
  vector<int> res;
  const vector<int> nchild = nchildvec();
  for (ENode n = 0; n < nodes(); ++n)
    if (nchild[n] > 0)
      res.push_back (n);
  return res;
}

vector<ENode> ETree::children (ENode node) const
{
  vector<ENode> res;
  for (ENode n = node + 1; n < nodes(); ++n)
    if (parent[n] == node)
      res.push_back (n);
    else if (parent[n] < node)
      break;
  return res;
}

vector<ENode> ETree::ancestors (ENode node) const
{
  vector<ENode> res;
  for (ENode n = parent[node]; n >= 0; n = parent[n])
    res.push_back (n);
  vector<ENode> rev (res.rbegin(), res.rend());  // reverse the order, inelegantly
  return rev;
}

vector<ENode> ETree::strictDescendants (ENode node) const
{
  const ENode endnode = firstEclipsed (node);
  vector<ENode> res (endnode - node - 1);
  for (ENode n = node + 1; n < endnode; ++n)
    res[n - (node + 1)] = n;
  return res;
}

vector<ENode> ETree::descendants (ENode node) const
{
  const ENode endnode = firstEclipsed (node);
  vector<ENode> res (endnode - node);
  for (ENode n = node; n < endnode; ++n)
    res[n - node] = n;
  return res;
}

vector<ENode> ETree::nondescendants (ENode node) const
{
  vector<ENode> res;
  vector<bool> from (nodes(), FALSE);
  for (ENode n = 0; n < nodes(); ++n)
    if ((parent[n] >= 0 && from[parent[n]]) || n == node)
      from[n] = TRUE;
    else
      res.push_back (n);
  return res;
}

vector<ENode> ETree::eclipsed (ENode node) const
{
  vector<ENode> res;
  vector<bool> from (nodes(), FALSE);
  from[node] = TRUE;
  for (ENode n = node + 1; n < nodes(); ++n)
    if (from[parent[n]])
      from[n] = TRUE;
    else
      res.push_back (n);
  return res;
}

vector<ENode> ETree::all() const
{
  vector<ENode> res (nodes());
  for (ENode n = 0; n < nodes(); ++n)
    res[n] = n;
  return res;
}

vector<int> ETree::nchildvec() const
{
  vector<int> nchild (nodes(), (int) 0);
  for (ENode node = 1; node < nodes(); ++node)
    ++nchild[parent[node]];
  return nchild;
}

ETree ETree::subtree (const vector<int>& nodes) const
{
  ETree etree (0);
  map<int,int> new_node_index;  // map from old node index to new node index
  for_const_contents (vector<int>, nodes, n)
    {
      new_node_index[*n] = etree.parent.size();
      const int p = parent[*n];
      if (new_node_index.find (p) == new_node_index.end())
	etree.parent.push_back (-1);
      else
	etree.parent.push_back (new_node_index[p]);
    }
  return etree;
}

// ESubtrees

ESubtrees::ESubtrees (const ETree& tree)
  : leaves (tree.leaves()),
    internals (tree.internals()),
    all (tree.all()),
    firstEclipsed (tree.nodes()),
    children (tree.nodes()),
    ancestors (tree.nodes()),
    strictDescendants (tree.nodes()),
    descendants (tree.nodes()),
    nondescendants (tree.nodes()),
    eclipsed (tree.nodes())
{
  for (ENode n = 0; n < tree.nodes(); ++n)
    {
      firstEclipsed[n] = tree.firstEclipsed (n);
      children[n] = tree.children (n);
      ancestors[n] = tree.ancestors (n);
      strictDescendants[n] = tree.strictDescendants (n);
      descendants[n] = tree.descendants (n);
      nondescendants[n] = tree.nondescendants (n);
      eclipsed[n] = tree.eclipsed (n);
    }
}

// CollapsedETrans

basic_string<char> CollapsedETrans::etrans2string() const
{
  basic_string<char> s;
  if (null()) return s;
  s.append (src().estate2string());
  s.append ("=>");
  s.append (dest().estate2string());
  s.append (" (");
  for (int n = 0; n < nodes(); ++n)
    {
      if (n > 0)
	s.append (", ");
      s.append (tpath[n].tpath2string());
    }
  s.append (")");
  return s;
}

basic_string<char> CollapsedETrans::etrans2string (const TAlphabet& alphabet) const
{
  basic_string<char> s;
  if (null()) return s;
  s.append (src().estate2string (alphabet));
  s.append ("=>");
  s.append (dest().estate2string (alphabet));
  s.append (" (");
  for (int n = 0; n < nodes(); ++n)
    {
      if (n > 0)
	s.append (", ");
      s.append (tpath[n].tpath2string (alphabet));
    }
  s.append (")");
  return s;
}

// CollapsedESpaceSubsets

CollapsedESpaceSubsets::CollapsedESpaceSubsets (const TSpace& tspace, const ETree& tree)
{
  build (tspace, tree, 0, 0, FALSE);
}

void CollapsedESpaceSubsets::build (const TSpace& tspace, const ETree& tree,
				    ENode root, const TEmission* emission, bool wipeChildren)
{
  // check tree is valid
  if (!tree.valid())
    return;

  // count nodes & terminals
  const int nodes = tree.nodes();
  const int terminals = tspace.terminals;

  // initialise A, B
  A = vector<ABMap> (nodes);
  B = vector<ABMap> (nodes);

  // precompute ESubtrees
  const ESubtrees subtrees (tree);

  // precompute {MatchXY,InsertY,DeleteX} sets
  // should probably make a summary struct for these, like ESubtrees
  const vector<TState> insertY = tspace.allInsertY();
  const vector<TState> deleteX = tspace.allDeleteX();
  vector<TStateVec> matchXY (terminals);
  for (TTerm xterm = 0; xterm < terminals; ++xterm)
    matchXY[xterm] = tspace.allMatchX (xterm);

  // fill A and B in preorder traversal
  for (ENode node = subtrees.firstEclipsed[root] - 1; node >= root; --node)
    {
      // get children
      const vector<ENode>& child = subtrees.children[node];
      const int children = child.size();

      // Precompute a few set unions.
      // Amd[C][X] = union of {A[C][MatchXY]} and A[C][DeleteX]
      // Bmi[C][X] = union of {B[C][MatchXY]} and {B[C][InsertY]}
      // Bsi[C]    = union of B[C][Start] and {B[C][InsertY]}
      // AmdBmi[X] = union of Amd[C][X] and Bmi[C][X] for C = children - 1
      typedef vector<EStateList> EStateTermVec;
      vector<EStateTermVec> Amd (children, EStateTermVec (terminals));
      vector<EStateTermVec> Bmi (Amd);
      vector<EStateList> AmdBmi (terminals);
      vector<EStateList> Bsi (children);

      // loop through children
      for (int c = 0; c < children; ++c)
	{
	  // note if this is the last child
	  const bool lastchild = c == children - 1;

	  // make refs of A[child[c]] and B[child[c]]
	  ABMap& AC = A[child[c]];
	  ABMap& BC = B[child[c]];

	  // make ref of target Bsi
	  EStateList& BsiC = Bsi[c];
	  // get EStateList
	  const EStateList& BStart = BC[TransStart];
	  // add start to Bsi[C]
	  BsiC.insert (BsiC.end(), begin_end (BStart));

	  // loop through X-terminals
	  for (TTerm xterm = 0; xterm < terminals; ++xterm)
	    {
	      // make refs of target sets
	      EStateList& AmdCX = Amd[c][xterm];
	      EStateList& BmiCX = Bmi[c][xterm];
	      EStateList& AmdBmiX = AmdBmi[xterm];

	      // loop through Y-terminals, doing MatchXY
	      for (TTerm yterm = 0; yterm < terminals; ++yterm)
		{

		  // get EStateList's
		  const EStateList& AMatchXY = AC[matchXY[xterm][yterm]];
		  const EStateList& BMatchXY = BC[matchXY[xterm][yterm]];

		  // add matches to Amd[C][X] and Bmi[C][X]
		  AmdCX.insert (AmdCX.end(), begin_end (AMatchXY));
		  BmiCX.insert (BmiCX.end(), begin_end (BMatchXY));
		}

	      // add matches to AmdBmi[X]
	      if (lastchild)
		{
		  AmdBmiX.insert (AmdBmiX.end(), begin_end (AmdCX));
		  AmdBmiX.insert (AmdBmiX.end(), begin_end (BmiCX));
		}

	      // get EStateList's
	      const EStateList& ADeleteX = AC[deleteX[xterm]];
	      const EStateList& BInsertX = BC[insertY[xterm]];

	      // add insert/delete to Bsi[C] and Amd[C][X], and possibly delete to AmdBmi[X]
	      BsiC.insert (BsiC.end(), begin_end (BInsertX));
	      AmdCX.insert (AmdCX.end(), begin_end (ADeleteX));
	      if (lastchild)
		AmdBmiX.insert (AmdBmiX.end(), begin_end (ADeleteX));

	      // loop through Y-terminals, doing InsertY
	      for (TTerm yterm = 0; yterm < terminals; ++yterm)
		{
		  // get EStateList
		  const EStateList& BInsertY = BC[insertY[yterm]];
		  // add insert to Bmi[C][X] and possibly AmdBmi[X]
		  BmiCX.insert (BmiCX.end(), begin_end (BInsertY));
		  if (lastchild)
		    AmdBmiX.insert (AmdBmiX.end(), begin_end (BInsertY));
		}
	    }
	}

      // Initialise childStateList for combine()
      vector<EStateList*> childStateList (children, (EStateList*) 0);

      // Now come the calls to combine().
      // Fill A
      // TransMatchXY
      for (TTerm x = 0; x < terminals; ++x)
	{
	  for (int c = 0; c < children; ++c)
	    childStateList[c] = &Amd[c][x];
	  for (TTerm y = 0; y < terminals; ++y)
	    combine (tspace, tree, subtrees, emission, node, tspace.TransMatchXY(x,y), childStateList, A);
	}

      // TransWait, TransDeleteX
      for (int c = 0; c < children; ++c)
	childStateList[c] = &A[child[c]][TransWait];
      combine (tspace, tree, subtrees, emission, node, TransWait, childStateList, A);
      for (TTerm x = 0; x < terminals; ++x)
	combine (tspace, tree, subtrees, emission, node, tspace.TransDeleteX(x), childStateList, A);

      // TransStart
      for (int c = 0; c < children; ++c)
	childStateList[c] = &A[child[c]][TransStart];
      combine (tspace, tree, subtrees, emission, node, TransStart, childStateList, A);

      // TransEnd
      for (int c = 0; c < children; ++c)
	childStateList[c] = &A[child[c]][TransEnd];
      combine (tspace, tree, subtrees, emission, node, TransEnd, childStateList, A);

      // Fill B
      if (children == 0)   // handle leaf nodes separately
	for (TTerm x = 0; x < terminals; ++x)
	  combine (tspace, tree, subtrees, emission, node, tspace.TransInsertY(x), childStateList, B);  // TransInsertY
      else  // not a leaf node
	for (int emitchild = 0; emitchild < children; ++emitchild)
	  {
	    // Loop through X-terminals
	    for (TTerm x = 0; x < terminals; ++x)
	      {
		// set up match, insert
		childStateList[emitchild] = &Bmi[emitchild][x];
		for (int c = 0; c < emitchild; ++c)
		  childStateList[c] = &Amd[c][x];
		for (int c = emitchild + 1; c < children; ++c)
		  childStateList[c] = &A[child[c]][TransWait];
		
		// TransMatchXY
		for (TTerm y = 0; y < terminals; ++y)
		  combine (tspace, tree, subtrees, emission, node, tspace.TransMatchXY(x,y), childStateList, B);
		
		// TransInsertY
		if (emitchild == children - 1)
		  childStateList[emitchild] = &AmdBmi[x];
		combine (tspace, tree, subtrees, emission, node, tspace.TransInsertY(x), childStateList, B);
	      }

	      // set up start
	      childStateList[emitchild] = &Bsi[emitchild];
	      for (int c = 0; c < emitchild; ++c)
		childStateList[c] = &A[child[c]][TransStart];
	      for (int c = emitchild + 1; c < children; ++c)
		childStateList[c] = &A[child[c]][TransWait];
	      
	      // TransStart
	      combine (tspace, tree, subtrees, emission, node, tspace.TransStart, childStateList, B);
	  }
      
      // optionally wipe children, for (slight) space efficiency
      if (wipeChildren)
	for_const_contents (vector<ENode>, child, c)
	{
	  A[*c].clear();
	  B[*c].clear();
	}
    }
}

EStateList CollapsedESpaceSubsets::joint (const TSpace& tspace) const
{
  EStateList res;
  if (A.size())  // check that the tree has nodes
    {
      ABMap& Asubroot = (ABMap&) A[0];  // cast away const so we can use operator[]
      ABMap& Bsubroot = (ABMap&) B[0];  // cast away const so we can use operator[]
      res.insert (res.end(), begin_end (Asubroot[TransStart]));
      res.insert (res.end(), begin_end (Bsubroot[TransStart]));
      res.insert (res.end(), begin_end (Asubroot[TransEnd]));
      res.insert (res.end(), begin_end (Asubroot[TransWait]));
      for (TTerm y = 0; y < tspace.terminals; ++y)
	res.insert (res.end(), begin_end (Bsubroot[tspace.TransInsertY(y)]));
    }
  return res;
}

EStateList CollapsedESpaceSubsets::conditional (const TSpace& tspace) const
{
  EStateList res;
  if (A.size())  // check that the tree has nodes
    {
      ABMap& Asubroot = (ABMap&) A[0];  // cast away const so we can use operator[]
      ABMap& Bsubroot = (ABMap&) B[0];  // cast away const so we can use operator[]
      res.insert (res.end(), begin_end (Asubroot[TransStart]));
      res.insert (res.end(), begin_end (Bsubroot[TransStart]));
      res.insert (res.end(), begin_end (Asubroot[TransEnd]));
      res.insert (res.end(), begin_end (Asubroot[TransWait]));
      for (TTerm x = 0; x < tspace.terminals; ++x)
	for (TTerm y = 0; y < tspace.terminals; ++y)
	  {
	    res.insert (res.end(), begin_end (Asubroot[tspace.TransMatchXY(x,y)]));
	    res.insert (res.end(), begin_end (Bsubroot[tspace.TransMatchXY(x,y)]));
	  }
      for (TTerm x = 0; x < tspace.terminals; ++x)
	  res.insert (res.end(), begin_end (Asubroot[tspace.TransDeleteX(x)]));
      for (TTerm y = 0; y < tspace.terminals; ++y)
	  res.insert (res.end(), begin_end (Bsubroot[tspace.TransInsertY(y)]));
    }
  return res;
}

// EStateIterator

EStateIterator::~EStateIterator()
{ }

// EStateNullIterator

EStateNullIterator::EStateNullIterator()
{ }

bool EStateNullIterator::finished()
{
  return TRUE;
}

void EStateNullIterator::getCurrent (EState& estate) const
{
  return;
}

void EStateNullIterator::next()
{
  return;
}

void EStateNullIterator::reset()
{
  return;
}

// estateNullIterator

EStateNullIterator estateNullIterator;

// EStateNodeLabelIterator

EStateNodeLabelIterator::EStateNodeLabelIterator()
  : node (-1), label (TTermNull), done (TRUE)
{ }

EStateNodeLabelIterator::EStateNodeLabelIterator (ENode node, TState label)
  : node (node), label (label), done (FALSE)
{ }

bool EStateNodeLabelIterator::finished()
{
  return done;
}

void EStateNodeLabelIterator::getCurrent (EState& estate) const
{
  if (!done)
    estate.tstate[node] = label;
  return;
}

void EStateNodeLabelIterator::next()
{
  done = TRUE;
  return;
}

void EStateNodeLabelIterator::reset()
{
  done = FALSE;
  return;
}

// EStateUnionIterator

EStateUnionIterator::EStateUnionIterator()
  : setIter(), currentSetIter (setIter.begin())
{ }

EStateUnionIterator::EStateUnionIterator (const EStateIterVec& si)
  : setIter (si), currentSetIter (setIter.begin())
{
  skipEmpty();
}

bool EStateUnionIterator::finished()
{
  return currentSetIter == setIter.end();
}

void EStateUnionIterator::getCurrent (EState& estate) const
{
  if (currentSetIter != setIter.end())
    (*currentSetIter)->getCurrent (estate);
  return;
}

void EStateUnionIterator::next()
{
  if (currentSetIter != setIter.end())
    (*currentSetIter)->next();
  skipEmpty();
  return;
}

void EStateUnionIterator::reset()
{
  for_contents (EStateIterVec, setIter, si)
    (*si)->reset();
  currentSetIter = setIter.begin();
  skipEmpty();
  return;
}

void EStateUnionIterator::skipEmpty()
{
  while (currentSetIter != setIter.end() && (*currentSetIter)->finished())
    {
      (*currentSetIter)->reset();
      ++currentSetIter;
    }
}

// EStateProductIterator

EStateProductIterator::EStateProductIterator (const EStateIterVec& childIter)
  : childIter (childIter), done (FALSE)
{
  if (childIter.size() == 0)
    done = TRUE;
  else
    for_const_contents (EStateIterVec, childIter, ci)
      if ((*ci)->finished())
	{
	  done = TRUE;
	  break;
	}
}

bool EStateProductIterator::finished()
{
  return done;
}

void EStateProductIterator::getCurrent (EState& estate) const
{
  if (done)
    return;
  for_const_contents (EStateIterVec, childIter, ci)
    (*ci)->getCurrent (estate);
  return;
}

void EStateProductIterator::next()
{
  if (done)
    return;
  const int children = childIter.size();
  int c;
  for (c = 0; c < children; ++c)
    {
      childIter[c]->next();
      if (childIter[c]->finished())
	childIter[c]->reset();
      else
	break;
    }
  done = c == children;
  return;
}

void EStateProductIterator::reset()
{
  if (childIter.size())
    {
      for_contents (EStateIterVec, childIter, ci)
	(*ci)->reset();
      done = FALSE;
      for_const_contents (EStateIterVec, childIter, ci)
	if ((*ci)->finished())
	  {
	    done = TRUE;
	    break;
	  }
    }
  return;
}

// NodeLabelFilter

NodeLabelFilter::~NodeLabelFilter()
{ }

// PermissiveNodeLabelFilter

PermissiveNodeLabelFilter::PermissiveNodeLabelFilter()
{ }

bool PermissiveNodeLabelFilter::nodeLabelAllowed (ENode node, TState label) const
{
  return TRUE;
}

// ExclusiveNodeLabelFilter

ExclusiveNodeLabelFilter::ExclusiveNodeLabelFilter()
{ }

bool ExclusiveNodeLabelFilter::nodeLabelAllowed (ENode node, TState label) const
{
  return FALSE;
}

// permissiveNodeLabelFilter, exclusiveNodeLabelFilter

PermissiveNodeLabelFilter permissiveNodeLabelFilter;
ExclusiveNodeLabelFilter exclusiveNodeLabelFilter;

// EmissionNodeLabelFilter

EmissionNodeLabelFilter::EmissionNodeLabelFilter (const TSpace& tspace, const TEmission& emission)
  : tspace (tspace),
    emission (emission)
{ }

bool EmissionNodeLabelFilter::nodeLabelAllowed (ENode node, TState label) const
{
  return tspace.yterm (label) == emission.tterm[node];
}

// CollapsedEStateIterators::ENodeTStateIter

void CollapsedEStateIterators::ENodeTStateIter::makeAB()
{
  EStateIterVec Aprodvec, Bprodvec;
  for_contents (list<EStateProductIterator>, Aprod, p)
    Aprodvec.push_back (&*p);
  for_contents (list<EStateProductIterator>, Bprod, p)
    Bprodvec.push_back (&*p);
  A = EStateUnionIterator (Aprodvec);
  B = EStateUnionIterator (Bprodvec);
}

// CollapsedEStateIterators::ENodeIter

CollapsedEStateIterators::ENodeIter::ENodeIter (int terminals)
  : ttermIter (terminals)
{ }

int CollapsedEStateIterators::ENodeIter::terminals() const
{
  return ttermIter.size();
}

void CollapsedEStateIterators::ENodeIter::makeUnions (bool lastChild)
{
  // call makeAB()
  for_contents (TStateIterMap, tstateIter, ti)
    ti->second.makeAB();

  // create a temporary TSpace
  const TSpace tspace (terminals());
  
  // make a work variable for creating EStateUnionIterator's
  EStateIterVec itervec;

  // make Amd
  for (TTerm x = 0; x < terminals(); ++x)
    {
      itervec.clear();
      for (TTerm y = 0; y < terminals(); ++y)
	itervec.push_back (&tstateIter[tspace.TransMatchXY(x,y)].A);
      itervec.push_back (&tstateIter[tspace.TransDeleteX(x)].A);
      ttermIter[x].Amd = EStateUnionIterator (itervec);
    }

  // make Bmi
  for (TTerm x = 0; x < terminals(); ++x)
    {
      itervec.clear();
      for (TTerm y = 0; y < terminals(); ++y)
	{
	  itervec.push_back (&tstateIter[tspace.TransMatchXY(x,y)].B);
	  itervec.push_back (&tstateIter[tspace.TransInsertY(y)].B);
	}
      ttermIter[x].Bmi = EStateUnionIterator (itervec);
    }
 
  // make Bsi
  itervec.clear();
  itervec.push_back (&tstateIter[TransStart].B);
  for (TTerm y = 0; y < terminals(); ++y)
    itervec.push_back (&tstateIter[tspace.TransInsertY(y)].B);
  Bsi = EStateUnionIterator (itervec);

  // optionally make AmdBmi
  if (lastChild)
    for (TTerm x = 0; x < terminals(); ++x)
      {
	itervec.clear();
	itervec.push_back (&ttermIter[x].Amd);
	itervec.push_back (&ttermIter[x].Bmi);
	ttermIter[x].AmdBmi = EStateUnionIterator (itervec);
      }
}

// CollapsedEStateIterators

CollapsedEStateIterators::CollapsedEStateIterators (const TSpace& tspace, const ETree& tree, ENode root, const NodeLabelFilter& filter)
{
  build (tspace, tree, root, filter);
}

void CollapsedEStateIterators::build (const TSpace& tspace, const ETree& tree, ENode root, const NodeLabelFilter& filter)
{
  // check tree is valid
  if (!tree.valid())
    return;

  // count nodes & terminals
  const int nodes = tree.nodes();
  const int terminals = tspace.terminals;

  // initialise enodeIter
  enodeIter = vector<ENodeIter> (nodes, ENodeIter (terminals));

  // fill A and B in preorder traversal
  for (ENode node = tree.firstEclipsed (root) - 1; node >= root; --node)
    {
      // get children
      const vector<ENode> child = tree.children (node);
      const int children = child.size();

      // call makeUnions() on children
      if (children)
	{
	  for (int c = 0; c < children - 1; ++c)
	    enodeIter[child[c]].makeUnions (FALSE);
	  enodeIter[child.back()].makeUnions (TRUE);
	}

      // Initialise childStateIter and dummyChildStateIter for combine()
      EStateIterVec childStateIter (children, (EStateIterator*) 0);
      const EStateIterVec dummyChildStateIter;

      // Now come the calls to combine().
      // Fill A
      // TransMatchXY
      for (TTerm x = 0; x < terminals; ++x)
	{
	  for (int c = 0; c < children; ++c)
	    childStateIter[c] = &enodeIter[child[c]].ttermIter[x].Amd;
	  for (TTerm y = 0; y < terminals; ++y)
	    combine (tspace, tree, filter, node, tspace.TransMatchXY(x,y), childStateIter, FALSE);
	}

      // TransWait, TransDeleteX
      for (int c = 0; c < children; ++c)
	childStateIter[c] = &enodeIter[child[c]].tstateIter[TransWait].A;
      combine (tspace, tree, filter, node, TransWait, childStateIter, FALSE);
      for (TTerm x = 0; x < terminals; ++x)
	combine (tspace, tree, filter, node, tspace.TransDeleteX(x), childStateIter, FALSE);

      // TransStart
      for (int c = 0; c < children; ++c)
	childStateIter[c] = &enodeIter[child[c]].tstateIter[TransStart].A;
      combine (tspace, tree, filter, node, TransStart, childStateIter, FALSE);

      // TransEnd
      for (int c = 0; c < children; ++c)
	childStateIter[c] = &enodeIter[child[c]].tstateIter[TransEnd].A;
      combine (tspace, tree, filter, node, TransEnd, childStateIter, FALSE);

      // Fill B
      if (children == 0)   // handle leaf nodes separately
	{
	  // Loop through X-terminals
	  for (TTerm x = 0; x < terminals; ++x)
	    {
	      // TransInsertY
	      combine (tspace, tree, filter, node, tspace.TransInsertY(x), dummyChildStateIter, TRUE);
	      // TransMatchXY
	      for (TTerm y = 0; y < terminals; ++y)
		combine (tspace, tree, exclusiveNodeLabelFilter, node, tspace.TransMatchXY(x,y), dummyChildStateIter, TRUE);
	    }
	  // TransStart
	  combine (tspace, tree, exclusiveNodeLabelFilter, node, TransStart, dummyChildStateIter, TRUE);
	}
      else  // not a leaf node
	for (int emitchild = 0; emitchild < children; ++emitchild)
	  {
	    // Loop through X-terminals
	    for (TTerm x = 0; x < terminals; ++x)
	      {
		// set up match, insert
		childStateIter[emitchild] = &enodeIter[child[emitchild]].ttermIter[x].Bmi;
		for (int c = 0; c < emitchild; ++c)
		  childStateIter[c] = &enodeIter[child[c]].ttermIter[x].Amd;
		for (int c = emitchild + 1; c < children; ++c)
		  childStateIter[c] = &enodeIter[child[c]].tstateIter[TransWait].A;
		
		// TransMatchXY
		for (TTerm y = 0; y < terminals; ++y)
		  combine (tspace, tree, filter, node, tspace.TransMatchXY(x,y), childStateIter, TRUE);
		
		// TransInsertY
		if (emitchild == children - 1)
		  childStateIter[emitchild] = &enodeIter[child[emitchild]].ttermIter[x].AmdBmi;
		combine (tspace, tree, filter, node, tspace.TransInsertY(x), childStateIter, TRUE);
	      }

	      // set up start
	      childStateIter[emitchild] = &enodeIter[child[emitchild]].Bsi;
	      for (int c = 0; c < emitchild; ++c)
		childStateIter[c] = &enodeIter[child[c]].tstateIter[TransStart].A;
	      for (int c = emitchild + 1; c < children; ++c)
		childStateIter[c] = &enodeIter[child[c]].tstateIter[TransWait].A;

	      // TransStart
	      combine (tspace, tree, filter, node, TransStart, childStateIter, TRUE);
	  }
    }

  // Call makeUnions() for root
  // Treat the root as a last child, even though it may not be (if it's nonzero),
  // because this way all unions are conveniently set up for the root.
  enodeIter[root].makeUnions (TRUE);
}

EStateUnionIterator CollapsedEStateIterators::conditional() const
{
  EStateIterVec itervec;
  if (enodeIter.size())  // check that the tree has nodes
    {
      ENodeIter& subroot = (ENodeIter&) enodeIter[0];  // cast away const so we can use operator[]
      const TSpace tspace (subroot.terminals());  // create a temporary TSpace
      itervec.push_back (&subroot.tstateIter[TransStart].A);
      itervec.push_back (&subroot.tstateIter[TransEnd].A);
      itervec.push_back (&subroot.tstateIter[TransWait].A);
      itervec.push_back (&subroot.tstateIter[TransStart].B);
      for (TTerm x = 0; x < subroot.terminals(); ++x)
	{
	  itervec.push_back (&subroot.ttermIter[x].Amd);  // can't just use AmdBmi here, or B[InsertY] gets included multiple times
	  itervec.push_back (&subroot.tstateIter[tspace.TransInsertY(x)].B);
	  for (TTerm y = 0; y < subroot.terminals(); ++y)
	    itervec.push_back (&subroot.tstateIter[tspace.TransMatchXY(x,y)].B);
	}
    }
  return EStateUnionIterator (itervec);
}

EStateUnionIterator CollapsedEStateIterators::joint() const
{
  EStateIterVec itervec;
  if (enodeIter.size())  // check that the tree has nodes
    {
      ENodeIter& subroot = (ENodeIter&) enodeIter[0];  // cast away const so we can use operator[]
      itervec.push_back (&subroot.tstateIter[TransStart].A);
      itervec.push_back (&subroot.tstateIter[TransEnd].A);
      itervec.push_back (&subroot.tstateIter[TransWait].A);
      itervec.push_back (&subroot.Bsi);
    }
  return EStateUnionIterator (itervec);
}

// CollapsedEMatrix

CollapsedEMatrix::CollapsedEMatrix (const EStateList& estates, const TSpace& tspace, const ETree& tree)
{
  // create ESubtrees
  const ESubtrees subtrees (tree);
  // iterate through all possible transitions
  for_const_contents (EStateList, estates, src)
    for_const_contents (EStateList, estates, dest)
    {
      const CollapsedETrans current (*src, *dest, tspace, tree, subtrees);
      if (!current.null())
	etrans.insert (etrans.end(), current);
    }
}
