#ifndef EHMM_INCLUDED
#define EHMM_INCLUDED

#include <stdlib.h>
#include <set>
#include <map>
#include <list>
#include <string>
#include <vector>
#include "util/macros.h"

// Transducer states and terminals
typedef int TTerm;   // terminal
typedef int TState;  // full nonterminal (including {TransDeleteX,TransMatchXY} but not {TransDelete,TransMatch})

// A TAlphabet is a transducer terminal alphabet, i.e. a string of (case-redundant) alphabet symbols
typedef basic_string<char> TAlphabet;

// Vector of transducer states
typedef vector<TState> TStateVec;

// Transducer state space enumeration
struct TSpaceEnum
{
  // TTerm enumeration
  enum { TTermNull = -1 };  // used to indicate a null emission
  // TState enumeration
  enum { TransUndef = -1,
	 TransStart = 0, TransEnd = 1, TransWait = 2,
	 TransInsert = 3, TransMatch = 4, TransDelete = 5,
	 TransTypes = 6 };
  // Method to strip off XY info
  static inline TState type (TState state);  // strips off XY
  // Methods for testing if transitions between TState's are allowed
  static inline bool allowed (TState src, TState dest);
  // Conversion from TState<-->char
  static inline char tstate2char (TState state);
  static inline TState char2tstate (char c);
  // Conversion from TState to xy emit
  static inline unsigned int xy_emit (TState state);
};

// Transducer state space, with embedded terminals in {MatchXY,DeleteX,InsertY}
struct TSpace : TSpaceEnum
{
  // Terminal alphabet size
  int terminals;
  // Constructor
  TSpace (int terminals = 1);
  // TState creators
  inline TState TransInsertY (TTerm y) const;
  inline TState TransDeleteX (TTerm x) const;
  inline TState TransMatchXY (TTerm x, TTerm y) const;
  // internal lookup methods used by TState creators
  inline TState TransDeleteBegin() const;
  inline TState TransMatchBegin() const;
  inline TState TransMatchYBegin (int y = 0) const;
  inline TState TransTotalStates() const;
  // TState sets
  vector<TState> allInsertY() const;
  vector<TState> allDeleteX() const;
  vector<TState> allMatchXY() const;
  vector<TState> allMatchX (TTerm xterm) const;  // returns all MatchXY with X==xterm
  vector<TState> allMatchY (TTerm yterm) const;  // returns all MatchXY with Y==yterm
  // Breakdown of TState info
  inline TTerm xterm (TState state) const;  // gets X from {MatchXY,DeleteX}
  inline TTerm yterm (TState state) const;  // gets Y from {MatchXY,InsertY}
  // Conversion from TState-->string (appends yterm, or '.')
  inline basic_string<char> tstate2string (TState state, const TAlphabet& alphabet) const;
};

// Transducer path between two non-TransWait states, possibly via an intermediate TransWait-state
struct TPath : TSpaceEnum
{
  // data
  vector<TState> tstate;
  // constructors
  inline TPath();  // returns empty path
  inline TPath (TState src);  // returns single-element path, signifying no transition
  inline TPath (TState src, TState dest);  // returns empty path if no legal transition via TransWait
  // helpers
  inline bool null() const;  // TRUE if empty path
  // conversion to string
  basic_string<char> tpath2string() const;
  basic_string<char> tpath2string (const TAlphabet& alphabet) const;
};

// ENode's
typedef int ENode;

// Vector of ENode's
typedef vector<ENode> ENodeVec;

// ETree's
// A tree as defined in (Holmes, ISMB 2003) but with nodes offset by -2,
// The root node is -1, and the subroot is node 0. (The root node is ignored for most operations.)
// NB nodes must be sorted in preorder.
struct ETree
{
  // data
  vector<ENode> parent;  // The parent of each node, with parent[0] = -1
  // constructors
  ETree (int nodes = 1);  // initialises to a star phylogeny off the sub-root
  ETree (const ETree& etree);

  // Helpers
  inline int nodes() const;

  // Test methods
  // Methods for testing the whole tree
  bool valid() const;  // tests if this is a valid ETree (i.e. nodes sorted in depth-first preorder)
  bool binary() const;  // tests if this is a binary tree
  // Methods for testing node pairs
  inline bool ancestral (ENode anc, ENode desc) const;  // returns TRUE if anc is ancestral to desc

  // Method for finding first node eclipsed by a node
  inline ENode firstEclipsed (ENode node) const;

  // Subtrees and subforests
  // None of these include the root node, -1.
  // Sets relating to the entire tree.
  vector<ENode> leaves() const;  // leaf nodes
  vector<ENode> internals() const;  // internal nodes
  vector<ENode> all() const;  // returns all the nodes
  // Sets relating to individual nodes
  vector<ENode> children (ENode node) const;  // returns all the immediate children of node
  vector<ENode> ancestors (ENode node) const;  // returns all the ancestors of node (does not include the node itself)
  vector<ENode> strictDescendants (ENode node) const;  // returns all the descendants of node (does not include the node itself)
  vector<ENode> descendants (ENode node) const;  // returns all the descendants of node, including the node itself
  vector<ENode> nondescendants (ENode node) const;  // returns all nodes except the node & its descendants
  vector<ENode> eclipsed (ENode node) const;  // returns all the nodes eclipsed by node (see Holmes 2003)

  // Helper method for whole-tree operations
  vector<int> nchildvec() const;  // counts number of children of each node

  // method to create an ETree with the topology of a given subset of nodes
  ETree subtree (const vector<int>& nodes) const;

  // To/from string methods
  // etree2string() returns a New Hampshire tree (branches unlabeled) with optional leaf names A,B,C...
  //       (may generate non-ASCII or duplicate leaf names for large trees)
  // string2etree() discards node names and branch lengths, and is forgiving about the terminating semicolon.
  inline basic_string<char> etree2string (bool leafnames = FALSE) const;
  static inline ETree string2etree (const basic_string<char>& s);
};

// Summary structure for ETree subsets
struct ESubtrees
{
  // data
  vector<ENode> leaves, internals, all, firstEclipsed;
  vector<ENodeVec> children, ancestors, strictDescendants, descendants, nondescendants, eclipsed;  // indexed by ENode
  // constructors
  ESubtrees() { }
  ESubtrees (const ETree& tree);
  // accessor
  inline int nodes() const;
};

// Vector of possibly-null terminals (emission)
struct TEmission : TSpaceEnum
{
  // data
  vector<TTerm> tterm;
  // constructor
  inline TEmission (int nodes);
  // helpers
  inline int nodes() const;
  inline ENode emitter() const;
  // summarize(...) returns information pertinent to Felsenstein pruning implementations:
  //  deleters[]  ...these nodes need to be initialized with "delete" vectors
  //  matchers[]  ...these nodes (sorted postorder) propagate messages up the tree, multiplied by "match" matrices
  //  inserter    ...this node is the clique root; the Felsenstein likelihood is a product of incoming messages with the "insert" vector
  inline void summarize (const ETree& etree, ENode& inserter, vector<ENode>& deleters, vector<ENode>& matchers) const;
};


// EHMM states
struct EState : TSpaceEnum
{
  // data
  vector<TState> tstate;

  // constructors
  inline EState (int nodes = 0, TState init_tstate = TransUndef);
  inline EState (const EState& estate);  // copy constructor
  // relational operators
  inline bool operator== (const EState& estate) const;
  inline bool operator< (const EState& estate) const;

  // number of nodes
  inline int nodes() const;

  // special node accessors
  inline ENode active() const;  // finds the active node, i.e. max n for which tstate[n] != TransWait
  inline ENode emitter() const;  // finds the emitter, i.e. max n for which type(tstate[n]) == TransInsert

  // Start, end and wait states
  // Node subsets involved in detecting start and end states
  inline vector<ENode> endNodes() const;  // returns all n for which tstate[n] == TransEnd
  inline vector<ENode> startNodesFrom (ENode node) const;  // returns all n for which tstate[n] == TransStart && n >= node
  // allStates() method returns one of {TransState,TransEnd,TransWait}
  // if all tstate's are equal to one of those, otherwise TransUndef
  inline TState allStates() const;

  // Methods related to emit states
  // sync() method returns TRUE if this is an emit state,
  // i.e. emitter's descendants are xterm-yterm synchronised with parents
  inline bool sync (const TSpace& tspace, const ETree& tree) const;
  // Emissions & absorption
  inline TEmission emission (const TSpace& tspace) const;  // finds the emission, i.e. yterm's of emitter & descendants
  inline TTerm absorption (const TSpace& tspace) const;  // finds the absorption, i.e. xterm of root state (if emitter == 0) or null (if emitter > 0)
  // The emit set
  inline vector<ENode> emitSet() const;  // returns all n for which tstate[n] is one of {TransMatch,TransInsert} && n >= emitter
  inline vector<ENode> collapsedChanges() const;  // returns all nodes which change on entry to this state in the collapsed EHMM, i.e. all n for which tstate[n] is one of {TransMatch,TransInsert,TransDelete} && n >= emitter

  // Methods to convert to/from string
  inline basic_string<char> estate2string() const;
  static inline EState string2estate (const basic_string<char>& s);
  // Conversion to long string
  inline basic_string<char> estate2string (const TAlphabet& alphabet) const;
};

// List of EState's
typedef list<EState> EStateList;

// Collapsed EHMM transitions
struct CollapsedETrans : TSpaceEnum
{
  // Path data
  vector<TPath> tpath;
  // Constructor
  // Returns a null CollapsedETrans unless the transition is allowed in the collapsed EHMM state space.
  inline CollapsedETrans (const EState& src, const EState& dest, const TSpace& tspace,
			  const ETree& tree, const ESubtrees& subtrees);
  // helpers
  inline bool null() const;  // TRUE if tpath.size() == 0
  inline int nodes() const;
  inline EState src() const;
  inline EState dest() const;
  // conversion to string
  basic_string<char> etrans2string() const;
  basic_string<char> etrans2string (const TAlphabet& alphabet) const;
};

// Subsets used in construction of collapsed EHMM state space
// Once this is fully functional, I want to replace the explicit EStateList's with complex iterators.
struct CollapsedESpaceSubsets : TSpaceEnum
{
  // Typedefs
  typedef map<TState,EStateList> ABMap;

  // Data
  vector<ABMap> A, B;  // indexed by ENode and then TState, e.g. A[node][tstate]

  // Constructor
  CollapsedESpaceSubsets (const TSpace& tspace, const ETree& tree);

  // Build method.
  // Here C = *emission is an optional constraining TEmission. If emission != 0, then:
  //  (1) Let E = root.
  //  (2) Only descendant nodes N of E (including E itself) are visited, so only partial EState's are created.
  //  (3) At node N, the only partial EState's {P} that are created have P.emission(tspace)[N] = C[N].
  //        (this is enforced in the combine() method)
  //  (4) The set Q = { P: P.emission(tspace) == C } is thus given by B[E][tspace.TransInsertY(C[E])].
  // To find outgoing and incoming sets, we set E = C.emitter(). Then:
  //  For the outgoing set {P} of all destination EStates with TEmission C that can be entered from source EState S:
  //   Start with Q and set P[N] = S[N] for all N < E.
  //     (We need to incorporate this in combine(), so that illegal states are not created in this way.)
  //     (This means adding an optional constraining EState prefix,
  //		in additon to the optional constraining TEmission.)
  //     (We'll also want a rebuild() method, so we can do this for multiple EState prefixes.)
  //  For the incoming set {P} of all source EStates with TEmission C that can enter destination EState D:
  //   Let F = D.emitter(). Then either (i) F <= E; (ii) F is in C's absorb set; (iii) C is null; (iv) none of the above.
  //   If (i)-(iii), start with Q and set P[N] = D[N] for all N < E.
  //   If (iv), the set is empty.
  void build (const TSpace& tspace, const ETree& tree,
	      ENode root = 0, const TEmission* emission = 0, bool wipeChildren = FALSE);

  // Methods to return collapsed EHMM state spaces
  // Conditionally normalised EHMM allows all states at node zero
  // Jointly normalised EHMM allows only {TransStart,TransEnd,TransWait,TransInsert} at node zero
  EStateList conditional (const TSpace& tspace) const;
  EStateList joint (const TSpace& tspace) const;

  // helper method to get product of a vector of EStateLists
  // childStateList[c] gives EStateList for subtrees.descendants[c]
  inline static void combine (const TSpace& tspace,
			      const ETree& tree,
			      const ESubtrees& subtrees,
			      const TEmission* emission,
			      const ENode node,
			      const TState label,
			      const vector<EStateList*>& childStateList,
			      vector<ABMap>& abmap);
};

// EState iterators
// Virtual template
struct EStateIterator
{
  // pure virtual methods
  virtual bool finished() = 0;
  virtual void getCurrent (EState& estate) const = 0;  // writes partial EState label to estate
  virtual void next() = 0;
  virtual void reset() = 0;

  // virtual destructor
  virtual ~EStateIterator();
};

// Vector of pointers to EStateIterator's
typedef vector<EStateIterator*> EStateIterVec;

// Null iterator
struct EStateNullIterator : EStateIterator
{
  // constructor
  EStateNullIterator();

  // EStateIterator virtuals
  bool finished();
  void getCurrent (EState& estate) const;
  void next();
  void reset();
};

// Singleton instance of EStateNullIterator, used as a default by CollapsedEStateIterators
extern EStateNullIterator estateNullIterator;

// Node labels
class EStateNodeLabelIterator : public EStateIterator, public TSpaceEnum
{
private:
  // data
  ENode node;
  TState label;
  bool done;

public:
  // constructors
  EStateNodeLabelIterator();
  EStateNodeLabelIterator (ENode node, TState label);

  // EStateIterator virtuals
  bool finished();
  void getCurrent (EState& estate) const;
  void next();
  void reset();
};

// Unions
class EStateUnionIterator : public EStateIterator
{
private:
  // typedefs
  typedef EStateIterVec::const_iterator EStateSetSelector;

  // data
  EStateIterVec setIter;
  EStateSetSelector currentSetIter;

public:
  // constructors
  EStateUnionIterator();
  EStateUnionIterator (const EStateIterVec& setIter);

  // EStateIterator virtuals
  bool finished();
  void getCurrent (EState& estate) const;
  void next();
  void reset();

  // helper method to skip empty sets
  void skipEmpty();
};

// Subtree combinations
class EStateProductIterator : public EStateIterator
{
private:
  // data
  EStateIterVec childIter;
  bool done;

public:
  // constructor
  EStateProductIterator (const EStateIterVec& childIter);

  // EStateIterator virtuals
  bool finished();
  void getCurrent (EState& estate) const;
  void next();
  void reset();
};

// NodeLabelFilter is an abstract filter on node labels
struct NodeLabelFilter
{
  // pure virtual method
  virtual bool nodeLabelAllowed (ENode node, TState label) const = 0;
  // virtual destructor
  virtual ~NodeLabelFilter();
};

// PermissiveNodeLabelFilter is a NodeLabelFilter that lets anything through
struct PermissiveNodeLabelFilter : NodeLabelFilter
{
  // constructor
  PermissiveNodeLabelFilter();
  // NodeLabelFilter virtual method
  bool nodeLabelAllowed (ENode node, TState label) const;
};

// ExclusiveNodeLabelFilter is a NodeLabelFilter that doesn't let anything through
struct ExclusiveNodeLabelFilter : NodeLabelFilter
{
  // constructor
  ExclusiveNodeLabelFilter();
  // NodeLabelFilter virtual method
  bool nodeLabelAllowed (ENode node, TState label) const;
};

// Singleton instances of PermissiveNodeLabelFilter and ExclusiveNodeLabelFilter, used as defaults by CollapsedEStateIterators
extern PermissiveNodeLabelFilter permissiveNodeLabelFilter;
extern ExclusiveNodeLabelFilter exclusiveNodeLabelFilter;

// EmissionNodeLabelFilter is a NodeLabelFilter that only allows node labels consistent with a given TEmission
class EmissionNodeLabelFilter : public NodeLabelFilter
{
private:
  // data
  TSpace tspace;
  TEmission emission;

public:
  // constructor
  EmissionNodeLabelFilter (const TSpace& tspace, const TEmission& emission);
  // NodeLabelFilter virtual method
  bool nodeLabelAllowed (ENode node, TState label) const;
};


// CollapsedEStateIterators is a structure that holds all the info needed to make an iterator over collapsed EStates
// The EStateIterator's use a lot of internal pointers (see above prototypes), so we disallow copy constructors, etc.
class CollapsedEStateIterators : TSpaceEnum
{
private:
  // ENodeTStateIter structure holds all the EStateIterators for an (ENode,TState) pair
  struct ENodeTStateIter
  {
    // data
    EStateNodeLabelIterator nodeLabel;
    list<EStateProductIterator> Aprod, Bprod;
    EStateUnionIterator A, B;  // these point to Aprod, Bprod after a call to makeAB()
    // builder
    void makeAB();  // points A, B to Aprod, Bprod
  };

  // ENodeTTermIter structure holds all the EStateIterators for an (ENode,TTerm) pair
  struct ENodeTTermIter
  {
    // data
    EStateUnionIterator Amd, Bmi, AmdBmi;  // these point to ENodeTStateIter::{A,B} after a call to ENodeIter::makeUnions()
  };

  // ENodeIterators structure holds all the EStateIterators for an ENode
  typedef map<TState,ENodeTStateIter> TStateIterMap;
  struct ENodeIter
  {
    // data
    TStateIterMap tstateIter;  // indexed by TState
    vector<ENodeTTermIter> ttermIter;  // indexed by TTerm
    EStateUnionIterator Bsi;  // points to tstateIter[].B after a call to makeUnions()
    // constructor
    ENodeIter (int terminals = 0);
    // helper
    int terminals() const;
    // builder
    // this calls makeAB() for all tstateIter's, then sets up remaining unions, including AmdBmi if lastChild==TRUE
    // Note: this method is NOT careful about preserving the tstateIter map.
    // If it tries to access tstateIter[S] and no such entry exists, then a new entry will be created,
    // corrupting earlier iterators.
    // It is therefore essential that this function is not called until the subtrees are already built.
    void makeUnions (bool lastChild);
  };

  // private methods
  // build() method
  void build (const TSpace& tspace, const ETree& tree, ENode root, const NodeLabelFilter& filter);

  // combine() builder method
  inline void combine (const TSpace& tspace,
		       const ETree& tree,
		       const NodeLabelFilter& filter,
		       const ENode node,
		       const TState label,
		       const EStateIterVec& childStateIter,
		       bool descendantsInsert);  // products added to (descendantsInsert ? Bprod : Aprod)

  // disallow copy constructor by making it private
  CollapsedEStateIterators (const CollapsedEStateIterators& copy);

public:
  // data
  vector<ENodeIter> enodeIter;  // indexed by ENode

  // constructor
  CollapsedEStateIterators (const TSpace& tspace, const ETree& tree,
			    ENode root = 0,
			    const NodeLabelFilter& filter = permissiveNodeLabelFilter);

  // Methods to return iterators over collapsed EHMM state spaces.
  // Note: only one iterator should ever be created from a single CollapsedEStateIterators object!
  // Otherwise iterators will be invalidated.
  // Conditionally normalised EHMM allows all states at node zero
  // Jointly normalised EHMM allows only {TransStart,TransEnd,TransWait,TransInsert} at node zero
  EStateUnionIterator conditional() const;
  EStateUnionIterator joint() const;
};

// Transition matrix for collapsed EHMM
struct CollapsedEMatrix
{
  // data
  vector<CollapsedETrans> etrans;
  // constructor
  CollapsedEMatrix (const EStateList& estates, const TSpace& tspace, const ETree& tree);
};

#include "handel/inline_ehmm.h"

#endif /* EHMM_INCLUDED */
