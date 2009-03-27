#ifndef INLINE_EHMM_INCLUDED
#define INLINE_EHMM_INCLUDED

#include <stack>
#include "util/vector_output.h"

// Inline method defs for file "ehmm.h"

// TSpaceEnum

TState TSpaceEnum::type (TState state)
{
  return state % TransTypes;
}

unsigned int TSpaceEnum::xy_emit (TState state)
{
  switch (type (state))
    {
    case TransMatch: return 3;
    case TransInsert: return 2;
    case TransDelete: return 1;
    case TransWait: return 0;
    default: break;
    }
  // Should probably throw an exception here, but this module is independent of util/logfile.h
  // THROWEXPR ("Can't get emit profile for state " << (int) state);
  return 4;  // dummy value
}

bool TSpaceEnum::allowed (TState src, TState dest)
{
  // strip off XY info
  const TState srctype = type (src);
  const TState desttype = type (dest);
  // check compatibility
  switch (srctype)
    {
    case TransStart:
    case TransMatch:  
    case TransDelete: 
    case TransInsert:
      return desttype == TransInsert || desttype == TransWait;
    case TransWait:
      return desttype == TransMatch || desttype == TransDelete || desttype == TransEnd;
    case TransEnd:
    default:
      return FALSE;
    }
  return FALSE;  // unreachable
}

char TSpaceEnum::tstate2char (TState state)
{
  char c;
  switch (type (state))
    {
    case TransStart:  c = 'S'; break;
    case TransEnd:    c = 'E'; break;
    case TransInsert: c = 'I'; break;
    case TransWait:   c = 'W'; break;
    case TransMatch:  c = 'M'; break;
    case TransDelete: c = 'D'; break;
    default:          c = '?'; break;
    }
  return c;
}

TState TSpaceEnum::char2tstate (char c)
{
  TState state;
  switch (toupper(c))
    {
    case 'S': state = TransStart;  break;
    case 'E': state = TransEnd;    break;
    case 'I': state = TransInsert; break;
    case 'W': state = TransWait;   break;
    case 'M': state = TransMatch;  break;
    case 'D': state = TransDelete; break;
    default:  state = TransUndef;  break;
    }
  return state;
}

// TSpace

TState TSpace::TransInsertY (TTerm y) const
{
  return TransInsert + y * TransTypes;
}

TState TSpace::TransDeleteX (TTerm x) const
{
  return TransDelete + x * TransTypes * terminals;
}

TState TSpace::TransMatchXY (TTerm x, TTerm y) const
{
  return TransMatch + (y + x * terminals) * TransTypes;
}

TTerm TSpace::xterm (TState state) const
{
  TTerm t;
  switch (type (state))
    {
    case TransMatch:
    case TransDelete:
      t = state / (TransTypes * terminals);
      break;
    case TransStart:
    case TransWait:
    case TransInsert:
    case TransEnd:
    default:
      t = TTermNull;
      break;
    }
  return t;
}

TTerm TSpace::yterm (TState state) const
{
  TTerm t;
  switch (type (state))
    {
    case TransMatch:
    case TransInsert:
      t = (state / TransTypes) % terminals;
      break;
    case TransStart:
    case TransWait:
    case TransDelete:
    case TransEnd:
    default:
      t = TTermNull;
      break;
    }
  return t;
}

basic_string<char> TSpace::tstate2string (TState state, const TAlphabet& alphabet) const
{
  basic_string<char> s;
  const int A = alphabet.size();
  if (A != terminals)  // balk if the alphabet doesn't match our terminal size
    return s;
  const TState t = type (state);
  s.push_back (tstate2char (t));
  if (A > 1)
    {
      const TTerm x = xterm (state);
      const TTerm y = yterm (state);
      if (x != TTermNull)
	s.push_back (toupper (alphabet[x]));
      if (y != TTermNull)
	s.push_back (toupper (alphabet[y]));
    }
  return s;
}

// ETree

int ETree::nodes() const
{
  return parent.size();
}

bool ETree::ancestral (ENode anc, ENode desc) const
{
  for (ENode n = parent[desc]; n >= 0; n = parent[n])
    if (n == anc)
      return TRUE;
  return FALSE;
}

ENode ETree::firstEclipsed (ENode node) const
{
  for (ENode n = node + 1; n < nodes(); ++n)
    if (parent[n] < node)
      return n;
  return nodes();
}

basic_string<char> ETree::etree2string (bool leafnames) const
{
  basic_string<char> s;
  if (!valid()) { s = "INVALID"; return s; }  // must be a valid ETree
  stack<int> anc;
  char leafchar = 'A';
  for (ENode n = 0; n < nodes(); ++n)
    {
      const ENode p = parent[n];
      if (p == n - 1)  // first child?
	{
	  s.push_back ('(');
	  anc.push (p);
	}
      else  // previous node must have been a leaf
	{
	  if (leafnames)
	    s.push_back (leafchar++);
	  while (p < anc.top())  // close brackets
	    {
	      anc.pop();
	      s.push_back (')');
	    }
	  s.push_back (',');
	}
    }
  if (leafnames)
    s.push_back (leafchar++);  // last node is always a leaf
  while (anc.size())  // close all remaining brackets
    {
      anc.pop();
      s.push_back (')');
    }
  s.push_back (';');

  return s;
}

ETree ETree::string2etree (const basic_string<char>& s)
{
  ETree tree (0);
  stack<int> anc;
  for_const_contents (basic_string<char>, s, c)
    switch (*c)
      {
      case '(':  // opening bracket
	anc.push (tree.nodes() - 1);
	// fall through to comma case
      case ',':  // comma
	tree.parent.push_back (anc.top());
	break;
      case ')':  // closing bracket
	if (anc.empty())  // if too many closing brackets, return an invalid ETree
	  return ETree (0);
	anc.pop();
	break;
      case ';':  // terminating semicolon
	if (!anc.empty())  // if too few closing brackets, return an invalid ETree
	  return ETree (0);
	return tree;  // ignore rest of string and return
	break;
      default:  // ignore other characters
	break;
      }
  // No ';' terminator found, but return anyway
  if (!anc.empty())  // if too few closing brackets, return an invalid ETree
    return ETree (0);
  return tree;
}

// ESubtrees

int ESubtrees::nodes() const
{
  return all.size();
}

// TEmission

TEmission::TEmission (int nodes)
  : tterm (nodes, (TTerm) TTermNull)
{ }

int TEmission::nodes() const
{
  return tterm.size();
}

ENode TEmission::emitter() const
{
  for (ENode n = 0; n < nodes(); ++n)
    if (tterm[n] != TTermNull)
      return n;
  return -1;
}

void TEmission::summarize (const ETree& etree,
			   ENode& inserter,
			   vector<ENode>& deleters,
			   vector<ENode>& matchers) const
{
  inserter = emitter();
  deleters.clear();
  matchers.clear();
  for (ENode n = etree.nodes() - 1; n > inserter; --n)
    {
      const ENode p = etree.parent[n];
      const TTerm ttn = tterm[n];
      const TTerm ttp = p >= 0 ? tterm[p] : TTermNull;
      if (ttp != TTermNull)
	if (ttn == TTermNull)
	  deleters.push_back (n);
	else
	  matchers.push_back (n);
    }
}

// EState

EState::EState (int nodes, TState init_tstate)
  : tstate (nodes, init_tstate)
{ }

EState::EState (const EState& estate)
  : tstate (estate.tstate)
{ }

bool EState::operator== (const EState& estate) const
{
  return tstate == estate.tstate;
}

bool EState::operator< (const EState& estate) const
{
  return tstate < estate.tstate;
}

int EState::nodes() const
{
  return tstate.size();
}

ENode EState::active() const
{
  for (ENode n = nodes() - 1; n > 0; --n)
    if (tstate[n] != TransWait)
      return n;
  return 0;
}

ENode EState::emitter() const
{
  for (ENode n = nodes() - 1; n > 0; --n)
    if (type (tstate[n]) == TransInsert)
      return n;
  return 0;
}

vector<ENode> EState::endNodes() const
{
  vector<ENode> res;
  for (ENode n = 0; n < nodes(); ++n)
    if (tstate[n] == TransEnd)
      res.push_back (n);
  return res;
}

vector<ENode> EState::startNodesFrom (ENode node) const
{
  vector<ENode> res;
  for (ENode n = node; n < nodes(); ++n)
    if (tstate[n] == TransStart)
      res.push_back (n);
  return res;
}

TState EState::allStates() const
{
  if (tstate.size() == 0)
    return TransUndef;
  const TState s = tstate[0];
  if (s < TransStart || s > TransWait)
    return TransUndef;
  for (int n = 1; n < nodes(); ++n)
    if (tstate[n] != s)
      return TransUndef;
  return s;
}

bool EState::sync (const TSpace& tspace, const ETree& tree) const
{
  const ENode e = emitter();
  vector<TTerm> yterm (tree.nodes(), TTermNull);
  yterm[e] = tspace.yterm (tstate[e]);
  for (ENode n = e + 1; n < tree.nodes(); ++n)
    {
      const TTerm parenty = yterm[tree.parent[n]];
      const TState s = tstate[n];
      if (tspace.xterm(s) != parenty)
	return FALSE;
      yterm[n] = tspace.yterm(s);
    }
  return TRUE;
}

TEmission EState::emission (const TSpace& tspace) const
{
  TEmission res (tstate.size());
  for (ENode n = emitter(); n < nodes(); ++n)
    res.tterm[n] = tspace.yterm (tstate[n]);
  return res;
}

TTerm EState::absorption (const TSpace& tspace) const
{
  return emitter() > 0 || tstate.size() == 0 ? TTermNull : tspace.xterm (tstate[0]);
}

vector<ENode> EState::emitSet() const
{
  vector<ENode> res;
  for (ENode n = emitter(); n < nodes(); ++n)
    if (tstate[n] == TransMatch || tstate[n] == TransInsert)
      res.push_back (n);
  return res;
}

vector<ENode> EState::collapsedChanges() const
{
  // if any branches are in end state & we're collapsed, must be all end
  vector<ENode> ends = endNodes();
  if (ends.size())
    return ends;
  // return emitter & absorbing states
  vector<ENode> res;
  for (ENode n = emitter(); n < nodes(); ++n)
    if (tstate[n] == TransMatch
	|| tstate[n] == TransInsert
	|| tstate[n] == TransDelete)
      res.push_back (n);
  return res;
}

basic_string<char> EState::estate2string() const
{
  basic_string<char> s (nodes(), '?');
  const ENode e = emitter();
  for (ENode n = 0; n < e; ++n)
    s[n] = tolower (tstate2char (tstate[n]));
  for (ENode n = e; n < nodes(); ++n)
    {
      const TState ts = tstate[n];
      s[n] = tstate2char (ts);
      if (type(ts) != TransMatch && type(ts) != TransInsert)
	s[n] = tolower (s[n]);
    }
  return s;
}

EState EState::string2estate (const basic_string<char>& s)
{
  EState e (s.size());
  for (ENode n = 0; n < e.nodes(); ++n)
    e.tstate[n] = char2tstate (s[n]);
  return e;
}

basic_string<char> EState::estate2string (const TAlphabet& alphabet) const
{
  const int terminals = alphabet.size();
  const TSpace tspace (terminals);  // create a temporary TSpace
  basic_string<char> s;
  const TTerm rootterm = tspace.xterm (tstate[0]);
  const ENode e = emitter();
  if (terminals > 1)
    if (rootterm == TTermNull)
      s.push_back ('_');
    else
      if (e > 0)
	s.push_back (tolower (alphabet[rootterm]));
      else
	s.push_back (toupper (alphabet[rootterm]));
  for (ENode n = 0; n < e; ++n)
    {
      const TState ts = tstate[n];
      const TState tt = type (tstate[n]);
      const TState ty = tspace.yterm (ts);
      s.push_back (tolower (tstate2char (tt)));
      if (terminals > 1)
	if (tt == TransMatch || tt == TransInsert)
	  s.push_back (tolower (alphabet[ty]));
	else
	  s.push_back ('_');
    }
  for (ENode n = e; n < nodes(); ++n)
    {
      const TState ts = tstate[n];
      const TState tt = type (tstate[n]);
      const TState ty = tspace.yterm (ts);
      if (tt == TransMatch || tt == TransInsert)
	{
	  if (terminals > 1)
	    {
	      s.push_back (tolower (tstate2char (tt)));
	      s.push_back (toupper (alphabet[ty]));
	    }
	  else
	    s.push_back (toupper (tstate2char (tt)));
	}
      else
	{
	  s.push_back (tolower (tstate2char (tt)));
	  if (terminals > 1)
	    s.push_back ('_');
	}
    }
  return s;
}

// TPath

TPath::TPath()
  : tstate()
{ }

TPath::TPath (TState src)
  : tstate (1, src)
{ }

TPath::TPath (TState src, TState dest)
{
  if (allowed (src, dest))
    {
      tstate.push_back (src);
      tstate.push_back (dest);
    }
  else if (allowed (src, TransWait) && allowed (TransWait, dest))
    {
      tstate.push_back (src);
      tstate.push_back (TransWait);
      tstate.push_back (dest);
    }
}

bool TPath::null() const
{
  return tstate.size() == 0;
}

// CollapsedETrans

CollapsedETrans::CollapsedETrans (const EState& src, const EState& dest, const TSpace& tspace,
				  const ETree& tree, const ESubtrees& subtrees)
{
  // If this isn't a valid transition, we return without initialising tpath.
  // The formulae for transition probabilities are revised from those presented in
  // (Holmes, ISMB 2003), which are incorrect.
  // The correction is documented at the following URL:

  // http://biowiki.org/TransducerCompositionErrata


  // Get emitter of src & dest
  const ENode srcEmitter = src.emitter();
  const ENode destEmitter = dest.emitter();
  // Check that emitter is either nonincreasing, or is in the src's absorb or start set
  if (destEmitter > srcEmitter)
    {
      const TState destEmitterSrcState = src.tstate[destEmitter];
      if (destEmitterSrcState != TransStart && tspace.xterm (destEmitterSrcState) == TTermNull)
	return;
    }

  // Check that any subroot transition is strictly transducer-legal
  // (transitions at other nodes are allowed to use intermediate TransWait states)
  if (destEmitter == 0 && !allowed (src.tstate[0], dest.tstate[0]))
    return;

  // Find first nodes eclipsed by srcEmitter & destEmitter
  // If N >= srcEcl:  N is eclipsed by srcEmitter
  //    N >= destEcl: N is eclipsed by destEmitter
  //    N >= maxEcl:  N is eclipsed by both srcEmitter & destEmitter
  const ENode srcEcl = subtrees.eclipsed[srcEmitter].size() ? subtrees.eclipsed[srcEmitter].front() : tree.nodes();
  const ENode destEcl = subtrees.eclipsed[destEmitter].size() ? subtrees.eclipsed[destEmitter].front() : tree.nodes();
  const ENode maxEcl = srcEcl > destEcl ? srcEcl : destEcl;  // max(srcEcl,destEcl)

  // Check all nodes, n, for valid changes.
  // If 0 <= n < destEmitter, no change allowed.
  // If destEmitter <= n < destEcl, then n is descended from destEmitter and is allowed to change.
  // If destEcl <= n < srcEcl, then n is newly eclipsed (and must be TransWait in dest, if it's a valid EState).
  // If max(srcEcl,destEcl) <= n < nodes(), then n was eclipsed and stays eclipsed.

  // Only the first of these rules rejects all transitions:
  for (ENode n = 0; n < destEmitter; ++n)
    if (src.tstate[n] != dest.tstate[n])
      return;

  // Initialise new tpath (don't update yet, as some TPath's may still evaluate to null(), which means we should abort)
  vector<TPath> newtpath (tree.nodes());

  // Create the TPath's
  // First do the ones that can come out null()
  // destEmitter must change
  if ((newtpath[destEmitter] = TPath (src.tstate[destEmitter], dest.tstate[destEmitter])).null())
    return;
  
  // destEmitter < n < destEcl: change only if dest or src is nonwait
  for (ENode n = destEmitter + 1; n < destEcl; ++n)
    if (src.tstate[n] == TransWait && dest.tstate[n] == TransWait)
      newtpath[n] = TPath (TransWait);
    else if ((newtpath[n] = TPath (src.tstate[n], dest.tstate[n])).null())
      return;

  // destEcl <= n < srcEcl: change only if src is nonwait
  for (ENode n = destEcl; n < srcEcl; ++n)
    if (src.tstate[n] == TransWait)
      newtpath[n] = TPath (TransWait);
    else if ((newtpath[n] = TPath (src.tstate[n], dest.tstate[n])).null())
      return;

  // Next do the nodes that don't change
  // 0 <= n < destEmitter: no change
  for (ENode n = 0; n < destEmitter; ++n)
    newtpath[n] = TPath (src.tstate[n]);

  // max(srcEcl,destEcl) <= n < nodes(): no change
  for (ENode n = maxEcl; n < tree.nodes(); ++n)
    newtpath[n] = TPath (src.tstate[n]);
  
  // update tpath
  tpath.swap (newtpath);
}

bool CollapsedETrans::null() const
{
  return tpath.size() == 0;
}

int CollapsedETrans::nodes() const
{
  return tpath.size();
}

EState CollapsedETrans::src() const
{
  EState estate (tpath.size());
  for (int n = 0; n < estate.nodes(); ++n)
    {
      if (tpath[n].null())
	return EState();  // return a null EState, signfifying failure
      estate.tstate[n] = tpath[n].tstate.front();
    }
  return estate;
}

EState CollapsedETrans::dest() const
{
  EState estate (tpath.size());
  for (int n = 0; n < estate.nodes(); ++n)
    {
      if (tpath[n].null())
	return EState();  // return a null EState, signfifying failure
      estate.tstate[n] = tpath[n].tstate.back();
    }
  return estate;
}

// CollapsedESpaceSubsets

void CollapsedESpaceSubsets::combine (const TSpace& tspace,
				      const ETree& tree,
				      const ESubtrees& subtrees,
				      const TEmission* emission,
				      const ENode node,
				      const TState label,
				      const vector<EStateList*>& childStateList,
				      vector<ABMap>& abmap)
{
  // check the emission constraint
  if (emission)
    if (tspace.yterm (label) != emission->tterm[node])
      return;

  // figure out where to put the results
  EStateList& result = abmap[node][label];

  // create a work variable for new EStates, and set the new node label
  EState estate (tree.nodes());
  estate.tstate[node] = label;

  // count children
  const vector<ENode>& children = subtrees.children[node];
  const int nchild = children.size();
  // if no children (i.e. leaf), return a single EState
  if (nchild == 0)
    {
      result.push_back (estate);
      return;
    }

  // initialise iterators
  vector<EStateList::const_iterator> childBegin (nchild);
  vector<EStateList::const_iterator> childEnd (nchild);
  vector<EStateList::const_iterator> childIter (nchild);
  unsigned int newstates = 1;
  for (int c = 0; c < nchild; ++c)
    {
      childBegin[c] = childIter[c] = childStateList[c]->begin();
      childEnd[c] = childStateList[c]->end();
      newstates *= childStateList[c]->size();
    }
  // check if any child sets were zero
  if (newstates == 0)
    return;

  // loop through product of child spaces
  while (1)
    {
      // merge the child EState's
      for (int child = 0; child < nchild; ++child)
	{
	  const EState& childstate = *childIter[child];
	  for_const_contents (ENodeVec, subtrees.descendants[children[child]], desc)
	    estate.tstate[*desc] = childstate.tstate[*desc];
	}
      // store the new EState
      result.push_back (estate);

      // update the iterators
      int inc;
      for (inc = 0; inc < nchild; ++inc)
	if (++childIter[inc] == childEnd[inc])
	  childIter[inc] = childBegin[inc];
	else
	  break;
      if (inc == nchild)
	break;
    }

  // return
  return;
}

// CollapsedEStateIterators

void CollapsedEStateIterators::combine (const TSpace& tspace,
					const ETree& tree,
					const NodeLabelFilter& filter,
					const ENode node,
					const TState label,
					const EStateIterVec& childStateIter,
					bool descendantsInsert)
{
  // get the tstateIter map entry for this node label, or create it if it doesn't exist
  TStateIterMap& tstateIterMap (enodeIter[node].tstateIter);
  TStateIterMap::iterator tstateIterFind = tstateIterMap.find (label);
  if (tstateIterFind == tstateIterMap.end())
    {
      // annotate the node label for this entry
      tstateIterFind = tstateIterMap.insert (TStateIterMap::value_type (label, ENodeTStateIter())).first;
      tstateIterFind->second.nodeLabel = EStateNodeLabelIterator (node, label);
    }
  ENodeTStateIter& tstateIter (tstateIterFind->second);
  
  // create the EStateIterVec that we will actually use
  EStateIterVec prodStateIter;

  // check the filter
  if (filter.nodeLabelAllowed (node, label))
    {
      // add the node label itself to the start of prodStateIter
      prodStateIter.push_back (&tstateIter.nodeLabel);
      prodStateIter.insert (prodStateIter.end(), childStateIter.begin(), childStateIter.end());
    }
  else  // filter rejected, so stick a null iterator into prodStateIter
    prodStateIter.push_back (&estateNullIterator);

  // figure out whether we're using A or B
  list<EStateProductIterator>& prod (descendantsInsert ? tstateIter.Bprod : tstateIter.Aprod);

  // create a new product iterator and add it to A or B
  prod.insert (prod.end(), EStateProductIterator (prodStateIter));

  // return
  return;
}

#endif /* INLINE_EHMM_INCLUDED */
