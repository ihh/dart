#include "handel/transducer_sexpr.h"
#include "handel/transducer_keywords.h"
#include "handel/multiwaydp.h"
#include "seq/pkeywords.h"

using std::fixed;

// regexp for detecting algebraic operator keywords
Regexp tsexpr_keyword_re = "^(\\*|\\+|\\-|/|" TSEXPR_SUM ")$";

// regexp for detecting token-like strings
Regexp tsexpr_token_re = "^" TSEXPR_ESC_TOKEN_PREFIX "[0-9]+$";

// default singleton transducer state names
#define Singleton_start_state_name  "*S"
#define Singleton_end_state_name    "*E"
#define Singleton_insert_state_name "*I"
#define Singleton_wait_state_name   "*W"

// suffix that's appended to transducer names in order to disambiguate them
#define Transducer_name_disambiguator "+"


Transducer_SExpr_file::Transducer_SExpr_file (SExpr& transducer_file_sexpr)
  : alphabet_size (0), etree (0)
{
  // read alphabet
  SExpr* alph_sexpr = transducer_file_sexpr.find (TSEXPR_TOKEN);
  if (alph_sexpr)
    {
      for_const_contents (list<SExpr>, alph_sexpr->value().child, c)
	{
	  const sstring& sym = c->get_atom();
	  if (sym2tok.find (sym) != sym2tok.end())
	    THROWEXPR ("Duplicate alphabet symbol " << sym);
	  alphabet.push_back (sym);
	  sym2tok[sym] = alphabet_size++;
	}
    }

  // read transducers
  const vector<SExpr*> transducer_sexpr_vec = transducer_file_sexpr.find_all (TSEXPR_TRANSDUCER);
  map<sstring,sstring> transducer_renaming_map;
  for (int pass = 0; pass < 2; ++pass)  // do two passes; handle transducers named Singleton_transducer_name on second pass
    for_const_contents (vector<SExpr*>, transducer_sexpr_vec, transducer_sexpr_ptr)
      {
	// get transducer name
	SExpr& transducer_sexpr = **transducer_sexpr_ptr;
	const sstring& orig_transducer_name = transducer_sexpr(TSEXPR_NAME).get_atom();
	sstring transducer_name (orig_transducer_name);
	if (transducer_name == Singleton_transducer_name)
	  {
	    if (pass == 0)
	      continue;
	    do {
	      transducer_name << Transducer_name_disambiguator;
	    } while (trans_dict.find (transducer_name) != trans_dict.end());
	    CLOGERR << "WARNING: the special transducer name " << orig_transducer_name << " is reserved; changing to " << transducer_name << " to avoid ambiguity\n";
	  }
	else
	  if (pass > 0)
	    continue;
	transducer_renaming_map[orig_transducer_name] = transducer_name;

	if (trans_dict.find (transducer_name) != trans_dict.end())
	  THROWEXPR ("Duplicate transducer name " << transducer_name);

	CTAG(1,TRANSDUCER) << "Building transducer " << transducer_name << "\n";

	// get states
	sstring start_name, end_name;
	StateIndexMap state_index;
	map<sstring,int> state_type_map;
	vector<sstring> state_name, emit_label;
	vector<State_type> state_type;
	vector<array2d<PFunc> > pair_emit;

	const vector<SExpr*> state_sexpr_vec = transducer_sexpr.find_all (TSEXPR_STATE);
	for_const_contents (vector<SExpr*>, state_sexpr_vec, state_sexpr)
	  {
	    // get state name
	    const sstring& this_state_name = (**state_sexpr) (TSEXPR_NAME).get_atom();
	    if (state_index.find (this_state_name) != state_index.end())
	      THROWEXPR ("Duplicate state name " << this_state_name);

	    // figure out state type
	    const sstring& this_state_type_name = (**state_sexpr) (TSEXPR_TYPE).get_atom();
	    int this_state_type = (int) UndefinedState;
	    int this_state_index = (int) state_name.size();
	    int this_dim = 0;  // number of symbols we have to deal with
	    array2d<PFunc> this_pair_emit;

	    if (this_state_type_name == sstring(TSEXPR_START))
	      {
		if (start_name.size())
		  THROWEXPR ("Duplicate start state");
		start_name = this_state_name;
		this_state_type = this_state_index = (int) Start;
	      }

	    else if (this_state_type_name == sstring(TSEXPR_END))
	      {
		if (end_name.size())
		  THROWEXPR ("Duplicate end state");
		end_name = this_state_name;
		this_state_type = this_state_index = (int) End;
	      }

	    else if (this_state_type_name == sstring(TSEXPR_WAIT))
	      this_state_type = (int) TransducerWaitType;

	    else if (this_state_type_name == sstring(TSEXPR_INSERT))
	      {
		this_state_type = (int) TransducerInsertType;
		this_dim = 1;
	      }

	    else if (this_state_type_name == sstring(TSEXPR_DELETE))
	      {
		this_state_type = (int) TransducerDeleteType;
		this_dim = 1;
	      }

	    else if (this_state_type_name == sstring(TSEXPR_MATCH))
	      {
		this_state_type = (int) TransducerMatchType;
		this_dim = 2;
	      }

	    // figure out state label
	    if (this_dim)
	      this_pair_emit.resize (alphabet_size, alphabet_size, PFunc (1.));
	    SExpr* this_emit_label_sexpr = (**state_sexpr).find (TSEXPR_LABEL);
	    sstring this_emit_label;
	    if (this_emit_label_sexpr)
	      {
		if (!this_emit_label_sexpr->value().is_atom())
		  THROWEXPR ("Emit labels of branch transducers must be atoms, not compound expressions");
		this_emit_label = this_emit_label_sexpr->value().get_atom();
		if (this_dim == 0)
		  {
		    CTAG(8,TRANSDUCER_WARNINGS) << "Warning: ignoring emit label '" << this_emit_label << "' on non-emit state\n";
		    this_emit_label.clear();
		  }
		else
		  {
		    if (tsexpr_keyword_re.Match (this_emit_label.c_str()))
		      THROWEXPR ("Uh-uh. You can't use '" << this_emit_label << "' as an emit label: it's a reserved algebraic operator. Could cause confusion. Sorry.");
		    if (tsexpr_token_re.Match (this_emit_label.c_str()))
		      THROWEXPR ("Emit label '" << this_emit_label << "' looks a bit too much like a token for my comfort level. Can you use something else, please?");

		    Alphabet_group pgroup;
		    if (sym2pgroup.find (this_emit_label) == sym2pgroup.end())
		      {
			// create new PGroup
			pgroup = pscores.new_alphabet_group (alphabet, this_dim, this_emit_label.c_str(), false);
			sym2pgroup[this_emit_label] = pgroup;
			// munge group_suffix. this is a pretty heinous abuse of an already fucked-up member variable
			munge_group_suffix (pscores, pgroup, this_emit_label);
		      }
		    else
		      {
			pgroup = sym2pgroup[this_emit_label];
			// type-check against existing PGroup
			const int old_dim = sym2pgroup[this_emit_label].word_len;
			if (old_dim != this_dim)
			  THROWEXPR ("Dimension mismatch for repeated emit label '" << this_emit_label << "'");
		      }

		    // set up pair_emit
		    vector<int> sym (this_dim, (int) 0);
		    switch (this_state_type)
		      {
		      case TransducerMatchType:
			for (sym[0] = 0; sym[0] < alphabet_size; ++sym[0])
			  for (sym[1] = 0; sym[1] < alphabet_size; ++sym[1])
			    Pair_trans_funcs::static_match_val (this_pair_emit, sym[0], sym[1]) = pgroup[sym];
			break;

		      case TransducerInsertType:
			for (sym[0] = 0; sym[0] < alphabet_size; ++sym[0])
			  Pair_trans_funcs::static_insert_val (this_pair_emit, sym[0]) = pgroup[sym];
			break;

		      case TransducerDeleteType:
			for (sym[0] = 0; sym[0] < alphabet_size; ++sym[0])
			  Pair_trans_funcs::static_delete_val (this_pair_emit, sym[0]) = pgroup[sym];
			break;

		      default:
			THROWEXPR ("Encountered unexpected emit label");
			break;
		      }
		  }
	      }

	    // store
	    state_index[this_state_name] = this_state_index;
	    if (this_state_index >= 0)
	      {
		state_name.push_back (this_state_name);
		state_type.push_back ((State_type) this_state_type);
		emit_label.push_back (this_emit_label);
		pair_emit.push_back (this_pair_emit);
	      }
	    state_type_map[this_state_name] = this_state_type;
	  }

	// check for start & end states
	if (start_name.size() == 0)
	  THROWEXPR ("No start state");

	if (end_name.size() == 0)
	  THROWEXPR ("No end state");

	// add extra wait states if there are any direct transitions to non-insert states
	const vector<SExpr*> transition_sexpr_vec = transducer_sexpr.find_all (TSEXPR_TRANSITION);
	typedef map<int,int> DummyStateMap;
	DummyStateMap dummy_wait_state;
	for_const_contents (vector<SExpr*>, transition_sexpr_vec, transition_sexpr)
	  {
	    // get source, dest states
	    const sstring& src_state_name = (**transition_sexpr) (TSEXPR_FROM).get_atom();
	    if (state_index.find (src_state_name) == state_index.end())
	      THROWEXPR ("State not found: " << src_state_name);

	    const sstring& dest_state_name = (**transition_sexpr) (TSEXPR_TO).get_atom();
	    if (state_index.find (dest_state_name) == state_index.end())
	      THROWEXPR ("State not found: " << dest_state_name);

	    const int src_state = state_index[src_state_name];

	    const int src_type = state_type_map[src_state_name];
	    const int dest_type = state_type_map[dest_state_name];

	    // check for transitions from non-wait states to absorbing/end states
	    if (src_type != TransducerWaitType
		&& (dest_type == TransducerMatchType || dest_type == TransducerDeleteType || dest_type == TransducerEndType)
		&& dummy_wait_state.find (src_state) == dummy_wait_state.end())
	      {
		// create name & index for new wait state, and add to list of states
		const int new_wait_state_index = state_name.size();
		sstring new_wait_state_name;
		for (int n = 0; true; ++n)
		  {
		    new_wait_state_name = src_state_name;
		    new_wait_state_name << TSEXPR_WAIT_SUFFIX;
		    if (n > 0)
		      new_wait_state_name << n;
		    if (state_type_map.find (new_wait_state_name) == state_type_map.end())
		      break;
		  }

		CTAG(8,TRANSDUCER_WARNINGS) << "Adding 'dummy' wait state " << new_wait_state_name << " splitting transition from state " << src_state_name << " to state " << dest_state_name << '\n';

		state_name.push_back (new_wait_state_name);
		state_type.push_back ((State_type) TransducerWaitType);
		emit_label.push_back (sstring());
		pair_emit.push_back (array2d<PFunc>());
		state_type_map[new_wait_state_name] = TransducerWaitType;
		dummy_wait_state[src_state] = new_wait_state_index;
	      }
	  }

	// build transducer
	trans_dict[transducer_name] = Pair_transducer_funcs (state_name.size());
	Pair_trans_funcs& transducer = trans_dict[transducer_name];

	transducer.start_name = start_name;
	transducer.end_name = end_name;
	transducer.state_name = state_name;
	transducer.state_type = state_type;
	transducer.emit_label = emit_label;
	transducer.pair_emit = pair_emit;
	transducer.alphabet_size = alphabet_size;

	// get transitions
	for_const_contents (vector<SExpr*>, transition_sexpr_vec, transition_sexpr)
	  {
	    // get source, dest states
	    const sstring& src_state_name = (**transition_sexpr) (TSEXPR_FROM).get_atom();
	    const sstring& dest_state_name = (**transition_sexpr) (TSEXPR_TO).get_atom();

	    int src_state = state_index[src_state_name];
	    int dest_state = state_index[dest_state_name];

	    int src_type = state_type_map[src_state_name];
	    int dest_type = state_type_map[dest_state_name];

	    // check dest is accessible from src
	    if (dest_type == TransducerStartType)
	      THROWEXPR ("Illegal transition from state " << src_state_name << " to start state " << dest_state_name);

	    switch (src_type)
	      {
	      case TransducerEndType:
		THROWEXPR ("Illegal transition from end state " << src_state_name);
		break;

	      case TransducerStartType:
	      case TransducerInsertType:
	      case TransducerMatchType:
	      case TransducerDeleteType:
		if (dest_type == TransducerMatchType || dest_type == TransducerDeleteType || dest_type == TransducerEndType)
		  {
		    src_state = dummy_wait_state[src_state];
		    src_type = TransducerWaitType;
		  }
		break;

	      case TransducerWaitType:
		if (dest_type == TransducerWaitType || dest_type == TransducerInsertType)
		  THROWEXPR ("Illegal transition from wait state " << src_state_name << " to state " << dest_state_name);
		break;

	      default:
		THROWEXPR ("bad transition type");
		break;
	      }

	    // check for duplicate transitions
	    if (!transducer.transition (src_state, dest_state).is_null())
	      THROWEXPR ("Duplicate transition from state " << src_state_name << " to state " << dest_state_name);

	    // make PFunc
	    PFunc f = 1.;
	    sstring trans_label ("1");
	    int n_pgroup = -1;

	    // get label, if it exists
	    SExpr* trans_label_sexpr = (**transition_sexpr).find (TSEXPR_LABEL);
	    if (trans_label_sexpr)
	      {
		if (!trans_label_sexpr->value().is_atom())
		  THROWEXPR ("Transition labels of branch transducers must be atoms, not compound expressions");
		trans_label = trans_label_sexpr->value().get_atom();
		if (sym2pgroup.find (trans_label) != sym2pgroup.end())
		  THROWEXPR ("Label '" << trans_label << "' is an emit label & can't be re-used as a transition label");
		if (sym2pvar.find (trans_label) == sym2pvar.end())
		  {
		    const PGroup pgroup = pscores.new_group ("dummy_group_name", trans_label.c_str());  // dummy single-PVar PGroup
		    n_pgroup = pgroup.group_idx;
		    f = (sym2pvar[trans_label] = pgroup[0]);
		  }
		else
		  {
		    const PVar pvar = sym2pvar[trans_label];
		    n_pgroup = pvar.group_idx;
		    f = pvar;
		  }
	      }

	    // add transition
	    transducer.transition (src_state, dest_state) = f;

	    // log
	    CTAG(1,TRANSDUCER) << "Added transition:\n";
	    CL << " names: " << src_state_name << " " << dest_state_name << "\n";
	    CL << " types: " << src_type << " " << dest_type << "\n";
	    CL << " indices: " << src_state << " " << dest_state << "\n";
	    CL << " label: " << trans_label << "\n";
	    CL << " pgroup: " << n_pgroup << "\n";
	    CL << " pfunc: ";
	    f.show (CL);
	    CL << " = ";
	    pfunc2stream (CL, pscores, f);
	    CL << "\n";
	  }

	// add dummy wait transitions
	for_const_contents (DummyStateMap, dummy_wait_state, dummy_wait)
	  {
	    transducer.transition (dummy_wait->first, dummy_wait->second) = PFunc (1.);
	    CTAG(1,TRANSDUCER) << "Added dummy wait transition from state "
			       << transducer.get_state_name (dummy_wait->first)
			       << " to state "
			       << transducer.get_state_name (dummy_wait->second)
			       << '\n';
	  }
      }

  // read tree
  SExpr& tree_sexpr = transducer_file_sexpr.find_or_die (TSEXPR_BRANCH);
  build_node (tree_sexpr, -1, transducer_renaming_map);
  CTAG(1,TRANSDUCER) << "ETree parent vector: " << etree.parent << "\n";
  CL << "ETree string representation: " << etree.etree2string() << "\n";

  // read Redelings-Suchard proposal schedule
  SExpr* redsuch_sexpr = transducer_file_sexpr.find (TSEXPR_SCHEDULE);
  if (redsuch_sexpr)
    {
      // build branch name->index map
      map<sstring,int> branch_name2index;
      for (int n = 0; n < etree.nodes(); ++n)
	branch_name2index[branch_name[n]] = n;
      // get proposal schedule
      const vector<SExpr*> branch_id_lists = redsuch_sexpr->values();
      for_const_contents (vector<SExpr*>, branch_id_lists, branch_id_list)
	{
	  vector<int> branch;
	  for_const_contents (list<SExpr>, (**branch_id_list).child, branch_id)
	    {
	      const sstring& bn = branch_id->get_atom();
	      if (branch_name2index.find (bn) == branch_name2index.end())
		THROWEXPR ("Branch " << bn << " not found");
	      const int bi = branch_name2index[bn];
	      if (path.find (bi) != path.end())
		THROWEXPR ("Can't resample a branch whose state path is constrained");
	      if (old_path.find (bi) == old_path.end())
		THROWEXPR ("Can't resample a branch unless a previous state path is specified");
	      branch.push_back (bi);
	    }
	  proposal_branches.push_back (branch);
	}
    }

  // set up cliques
  setup_cliques();

  // set composite name
  setup_composite_name();

  // read "value" and "bit-value" definitions
  read_pvars (transducer_file_sexpr, TSEXPR_BITVALUE, true);
  read_pvars (transducer_file_sexpr, TSEXPR_VALUE, false);

  vector<sstring> undef_vars;
  for (int g = 0; g < pscores.groups(); ++g)
    for (int v = 0; v < pscores.group_size(g); ++v)
      if (defined_pvars.find (PVar (g, v)) == defined_pvars.end())
	{
	  undef_vars.push_back (pscores.group_suffix[g][v]);
	  pscores[PVar(g,v)] = 0;
	}
  if (undef_vars.size() > 0 && (defined_pvars.size() > 0 || prof.size() > 0 || path.size() > 0))
    CTAG(8,TRANSDUCER_WARNINGS) << "Warning -- the following labels were undefined and will evaluate to one:\n" << undef_vars << '\n';
}

void Transducer_SExpr_file::munge_group_suffix (PScores& pscores, PGroup& pg, const sstring& emit_label)
{
  // munge group_suffix. this is a pretty heinous abuse of an already fucked-up member variable
  for_contents (vector<sstring>, pscores.group_suffix[pg.group_idx], gs)
    {
      sstring new_gs;
      new_gs << '(' << emit_label << ' ' << *gs << ')';
      *gs = new_gs;
    }
}

void Transducer_SExpr_file::setup_composite_name()
{
  vector<sstring> composite_name_string (etree.nodes());
  for (int n = etree.nodes() - 1; n >= 0; --n)
    {
      composite_name_string[n] << branch_trans_name[n];

      const vector<ENode> children = etree.children (n);
      if (children.size())
	{
	  for (int ci = 0; ci < (int) children.size(); ++ci)
	    composite_name_string[n] << (ci > 0 ? ' ' : '(') << composite_name_string[children[ci]];
	  composite_name_string[n] << ')';
	}
    }
  composite_name.clear();
  if (etree.nodes())
    composite_name << '(' << composite_name_string[0] << ')';
}

Transducer_SExpr_file::Transducer_SExpr_file (const vector<sstring>& alphabet,
					      const PScores& pscores,
					      const ETree& etree,
					      const vector<Pair_trans_funcs>& branch_trans,
					      const NodeProfileMap& prof,
					      const NodePathMap& path,
					      const NodePathMap& old_path)
  : alphabet (alphabet),
    pscores (pscores),
    etree (etree),
    branch_trans (branch_trans),
    prof (prof),
    path (path),
    old_path (old_path)
{
  init();
}

Transducer_SExpr_file::Transducer_SExpr_file (const Transducer_SExpr_file& trans)
  : alphabet (trans.alphabet),
    pscores (trans.pscores),
    etree (trans.etree),
    branch_name (trans.branch_name),
    tape_name (trans.tape_name),
    branch_trans (trans.branch_trans),
    prof (trans.prof),
    path (trans.path),
    old_path (trans.old_path),
    proposal_branches (trans.proposal_branches),
    band_coeff (trans.band_coeff)
{
  init();
}

void Transducer_SExpr_file::init()
{
  alphabet_size = alphabet.size();
  for (int tok = 0; tok < alphabet_size; ++tok)
    sym2tok[alphabet[tok]] = tok;

  for (int g = 0; g < pscores.groups(); ++g)
    for (int v = 0; v < pscores.group_size (g); ++v)
      defined_pvars.insert (PVar (g, v));

  setup_cliques();
  autoname_tree();
  autobuild_trans_dict();
  setup_composite_name();
}

void Transducer_SExpr_file::autobuild_trans_dict()
{
  trans_dict.clear();
  for (int b = 0; b < (int) branch_trans.size(); ++b)
    trans_dict[branch_trans_name[b]] = branch_trans[b];
}

void Transducer_SExpr_file::read_pvars (SExpr& transducer_file_sexpr, const char* keyword, bool is_bitscore)
{
  const vector<SExpr*> transducer_def_vec = transducer_file_sexpr.find_all (keyword);
  for_const_contents (vector<SExpr*>, transducer_def_vec, transducer_def_ptr)
    {
      const vector<SExpr*> sexvec = (*transducer_def_ptr)->values();
      for_const_contents (vector<SExpr*>, sexvec, tds_ptr)
	{
	  SExpr& transducer_def_sexpr = **tds_ptr;
	  if (transducer_def_sexpr.child.size() != 2)
	    THROWEXPR ("Malformed " << keyword << " expression");

	  SExpr& var_sexpr = transducer_def_sexpr[0];
	  const PVar pvar = resolve_pvar (var_sexpr);
	  if (pvar.group_idx >= 0)
	    {
	      if (defined_pvars.find (pvar) != defined_pvars.end())
		THROWEXPR ("Multiple definitions of " << var_sexpr);

	      Score sc;
	      const sstring dbl_string = transducer_def_sexpr[1].get_atom();
	      if (dbl_string == sstring (PK_INFINITE))
		{
		  if (!is_bitscore)
		    THROWEXPR ("Illegal infinite value in assignment " << transducer_def_sexpr);
		  sc = -InfinityScore;
		}
	      else if (is_bitscore)
		sc =  -Bits2Score(dbl_string.to_double());
	      else
		sc = Prob2Score (dbl_string.to_nonneg_double_strict());

	      set_pvar (pvar, sc);
	    }
	}
    }
}

void Transducer_SExpr_file::set_pvar (const PVar& pvar, Score sc)
{
  defined_pvars.insert (pvar);
  pscores[pvar] = sc;
}

PVar Transducer_SExpr_file::resolve_pvar (SExpr& var_sexpr)
{
  PVar pvar;
  if (var_sexpr.is_atom())
    {
      const sstring& sym = var_sexpr.get_atom();
      if (sym2pgroup.find (sym) != sym2pgroup.end())
	THROWEXPR ("Attempt to define array variable '" << sym << "' as if it were a scalar");
      if (sym2pvar.find (sym) == sym2pvar.end())
	CTAG(8,TRANSDUCER_WARNINGS) << "Warning: ignoring definition of unused scalar variable '" << sym << "'\n";
      else
	pvar = sym2pvar[sym];
    }
  else
    {
      const vector<sstring> sym = var_sexpr.atoms_to_strings();
      if (sym2pvar.find (sym[0]) != sym2pvar.end())
	{
	  if (sym.size() == 1)
	    THROWEXPR ("Scalar variable '" << sym << "' is enclosed in brackets, like a subscripted array variable; you probably didn't mean to do that & I'm going to be totally uptight about it");
	  THROWEXPR ("Attempt to subscript scalar variable '" << sym << "' as if it were an array");
	}
      if (sym.size() < 2)
	THROWEXPR ("Malformed array element: " << var_sexpr);
      if (sym2pgroup.find (sym[0]) == sym2pgroup.end())
	CTAG(8,TRANSDUCER_WARNINGS) << "Warning: ignoring definition of unknown array element '" << var_sexpr << "'\n";
      else
	{
	  const Alphabet_group alph_group = sym2pgroup[sym[0]];
	  const int dim = alph_group.word_len;
	  if ((int) sym.size() != dim + 1)
	    THROWEXPR ("Expected " << dim << " subscripts in array element " << var_sexpr);
	  vector<int> tok (dim);
	  for (int i = 0; i < dim; ++i)
	    {
	      if (sym2tok.find (sym[i+1]) == sym2tok.end())
		THROWEXPR ("Unknown alphabet token " << sym[i+1] << " in array subscript " << var_sexpr);
	      tok[i] = sym2tok[sym[i+1]];
	    }
	  pvar = alph_group[tok];
	}
    }
  return pvar;
}

int Transducer_SExpr_file::resolve_token (const sstring& tok)
{
  if (sym2tok.find (tok) == sym2tok.end())
    THROWEXPR ("Unknown alphabet token " << tok);
  return sym2tok[tok];
}

void Transducer_SExpr_file::build_node (SExpr& node, int parent, const map<sstring,sstring>& transducer_renaming_map)
{
  // get transducer (mandatory)
  const sstring& orig_transducer_name = node(TSEXPR_TRANSDUCER).get_atom();
  map<sstring,sstring>::const_iterator transducer_renaming_iter = transducer_renaming_map.find (orig_transducer_name);
  const sstring& transducer_name = transducer_renaming_iter == transducer_renaming_map.end() ? orig_transducer_name : transducer_renaming_iter->second;
  if (trans_dict.find (transducer_name) == trans_dict.end())
    THROWEXPR ("Transducer not found: " << transducer_name);

  // update ETree
  const int n_node = etree.nodes();
  etree.parent.push_back (parent);

  // make node name
  sstring n_name;
  if (node.find (TSEXPR_TO))
    n_name << node(TSEXPR_TO).get_atom();
  else
    n_name << auto_tape_name (n_node);
  tape_name[n_node] = n_name;

  // make parent node name here, for want of any place more intuitive to put this code
  if (node.find (TSEXPR_FROM))
    {
      if (tape_name.find (parent) != tape_name.end())
	if (tape_name[parent] != node(TSEXPR_FROM).get_atom())
	  THROWEXPR ("Input tape " << node(TSEXPR_FROM).get_atom() << " is already named " << tape_name[parent]);
      tape_name[parent] = node(TSEXPR_FROM).get_atom();
    }
  else if (parent < 0 && tape_name.find (parent) == tape_name.end())
    tape_name[parent] << TSEXPR_TOKEN_PREFIX << '0';
  
  // make branch name
  sstring this_branch_name;
  if (node.find (TSEXPR_NAME))
    this_branch_name = node(TSEXPR_NAME).get_atom();
  else
    this_branch_name << auto_branch_name (n_node);
  branch_name.push_back (this_branch_name);

  // log
  CTAG(1,TRANSDUCER) << "Branch " << this_branch_name
		     << " to tree node #" << n_node
		     << " (" << n_name
		     << ") with parent #" << parent
		     << " (" << tape_name[parent]
		     << ") has transducer " << transducer_name << "\n";

  // update branch->transducer map
  const Pair_trans_funcs& btrans = trans_dict[transducer_name];
  branch_trans_name.push_back (transducer_name);
  branch_trans.push_back (btrans);

  // state path specified?
  SExpr* path_sexpr = node.find (TSEXPR_STATE);
  SExpr* old_path_sexpr = node.find (TSEXPR_OLD_STATE);
  if (path_sexpr && old_path_sexpr)
    THROWEXPR ("Can't specify both '" << TSEXPR_STATE << "' and '" << TSEXPR_OLD_STATE << "'");
  if (path_sexpr || old_path_sexpr)
    {
      SExpr* this_path_sexpr = path_sexpr ? path_sexpr : old_path_sexpr;
      vector<int> this_path;
      StateIndexMap state_index = btrans.state_index();
      for_const_contents (list<SExpr>, this_path_sexpr->value().child, state_name_sexpr)
	{
	  const sstring& state_name = state_name_sexpr->get_atom();
	  if (state_index.find (state_name) == state_index.end())
	    THROWEXPR ("Can't resolve state " << state_name_sexpr->get_atom());
	  const int state = state_index[state_name];
	  if (state >= 0 && btrans.state_type[state])
	    this_path.push_back (state);
	}

      if (path_sexpr)
	path[n_node] = this_path;
      else
	old_path[n_node] = this_path;
    }

  // banding coefficient specified?
  SExpr* band_coeff_sexpr = node.find (TSEXPR_BAND_COEFF);
  if (band_coeff_sexpr)
    band_coeff[n_node] = band_coeff_sexpr->value().get_atom().to_double();

  // sequence, profile or bit-profile specified?
  bool seq_only = true, bit_matrix = false;
  SExpr* seq_sexpr = node.find (TSEXPR_SEQUENCE);
  if (!seq_sexpr)
    {
      seq_only = false;
      seq_sexpr = node.find (TSEXPR_PROFILE);
    }
  if (!seq_sexpr)
    {
      bit_matrix = true;
      seq_sexpr = node.find (TSEXPR_BITPROFILE);
    }

  // if sequence, profile or bit-profile specified, then read the bit-profile (or similar object)
  if (seq_sexpr)
    {
      SExpr& prof_sexpr = seq_sexpr->value();
      Score_profile prof_sc;
      for_contents (list<SExpr>, prof_sexpr.child, ssm_sexpr)
	{
	  Symbol_score_map ssm;
	  if (ssm_sexpr->is_atom())
	    ssm[resolve_token (ssm_sexpr->get_atom())] = 0;
	  else
	    {
	      if (seq_only)
		THROWEXPR ("Distributions over alphabet tokens are not allowed in "
			   << TSEXPR_SEQUENCE << " expressions; try " << TSEXPR_PROFILE << " or " << TSEXPR_BITPROFILE);
	      for_contents (list<SExpr>, ssm_sexpr->child, ss_sexpr)
		{
		  const double val = (*ss_sexpr)[1].get_atom().to_double();
		  const Score val_sc = bit_matrix ? -Bits2Score(val) : Prob2Score(val);
		  if (val_sc > -InfinityScore)
		    ssm[resolve_token ((*ss_sexpr)[0].get_atom())] = val_sc;
		}
	    }
	  if (ssm.empty())
	    THROWEXPR ("No valid characters at position " << (prof_sc.size() + 1) << " of " << seq_sexpr->tag());
	  prof_sc.push_back (ssm);
	}
      prof[n_node] = prof_sc;
    }

  // reconstruct this node?
  if (node.find (TSEXPR_RECONS))
    nodes_to_reconstruct.insert (n_node);

  // process nested branches recursively
  const vector<SExpr*> child_sexpr_vec = node.find_all (TSEXPR_BRANCH);
  for_const_contents (vector<SExpr*>, child_sexpr_vec, child_sexpr_ptr)
    build_node (**child_sexpr_ptr, n_node, transducer_renaming_map);
}

sstring Transducer_SExpr_file::auto_tape_name (int node)
{
  sstring n_name;
  n_name << TSEXPR_TOKEN_PREFIX << node + 1;
  return n_name;
}

sstring Transducer_SExpr_file::auto_branch_name (int node)
{
  sstring branch_name;
  branch_name << tape_name[etree.parent[node]] << TSEXPR_BRANCH_INFIX << tape_name[node];
  return branch_name;
}

void Transducer_SExpr_file::autoname_tree()
{
  for (int node = -1; node < etree.nodes(); ++node)
    {
      if (tape_name.find (node) == tape_name.end())
	tape_name[node] = auto_tape_name (node);

      if (node >= 0)
	{
	  if ((int) branch_name.size() <= node)
	    branch_name.push_back (auto_branch_name (node));
	  else if (branch_name[node].size() == 0)
	    branch_name[node] = auto_branch_name (node);
	}
    }

  for (int b = 0; b < (int) branch_trans.size(); ++b)
    {
      const sstring& btn = branch_name[b];
      if (b >= (int) branch_trans_name.size())
	branch_trans_name.push_back (btn);
      else if (branch_trans_name[b].size() == 0)
	branch_trans_name[b] = btn;
    }
}

void Transducer_SExpr_file::show_defs (ostream& out)
{
  // print alphabet
  if (alphabet_size)
    out << '(' << TSEXPR_TOKEN << " (" << alphabet << "))\n\n";

  // print label definitions
  if (defined_pvars.size())
    {
      out << '(' << TSEXPR_BITVALUE << '\n';
      for_const_contents (set<PVar>, defined_pvars, pvar)
	{
	  out << " (";
	  pvar->show (out, &pscores.group_suffix, true);
	  out << ' ' << score_sexpr (pscores[*pvar]) << ")\n";
	}
      out << ")\n";
    }

  // do we have any PVar definitions?
  const bool show_bitvalues = defined_pvars.size() > 0;

  // print branch transducers
  for_contents (SymTrans, trans_dict, bt)
    {
      out << '\n';
      bt->second.show_sexpr (out, bt->first, &pscores, show_bitvalues);
    }
}

void Transducer_SExpr_file::show_tree (ostream& out, int base_indent, bool show_composite_paths)
{
  // print tree
  const sstring base_indent_string (base_indent, ' ');
  vector<sstring> branch_string (etree.nodes());
  for (int n = etree.nodes() - 1; n >= 0; --n)
    {
      int nested_indent = base_indent;
      for (int a = n; a > 0; a = etree.parent[a])
	++nested_indent;
      const sstring nested_indent_string (nested_indent, ' ');

      branch_string[n] << nested_indent_string
		       << '(' << TSEXPR_BRANCH
		       << " (" << TSEXPR_FROM << ' ' << tape_name[etree.parent[n]]
		       << ") (" << TSEXPR_TO << ' ' << tape_name[n]
		       << ") (" << TSEXPR_NAME << ' ' << branch_name[n]
		       << ") (" << TSEXPR_TRANSDUCER << ' ' << branch_trans_name[n]
		       << ')';

      bool show_this_path = true, resample_this_path = false;
      NodePathMap::const_iterator path_iter = path.find (n);
      if (path_iter == path.end())
	{
	  resample_this_path = true;
	  path_iter = old_path.find (n);
	  if (path_iter == old_path.end())
	    show_this_path = false;
	}

      if (show_this_path)
	{
	  branch_string[n] << '\n' << nested_indent_string << ' ';
	  show_path (branch_trans[n], path_iter->second, resample_this_path ? TSEXPR_OLD_STATE : TSEXPR_STATE, branch_string[n]);

	  branch_string[n] << '\n' << nested_indent_string << ' ';
	  show_path_types (branch_trans[n], path_iter->second, branch_string[n]);
	}

      if (band_coeff.find (n) != band_coeff.end())
	branch_string[n] << '\n' << nested_indent_string << " (" << TSEXPR_BAND_COEFF << ' ' << band_coeff[n] << ')';

      if (prof.find (n) != prof.end())
	{
	  branch_string[n] << '\n' << nested_indent_string << ' ';
	  show_prof (prof[n], branch_string[n], nested_indent + 1);
	}

      if (nodes_to_reconstruct.find (n) != nodes_to_reconstruct.end())
	branch_string[n] << '\n' << nested_indent_string << " (" << TSEXPR_RECONS << ')';

      const vector<ENode> children = etree.children (n);
      if (children.size())
	{
	  for (int ci = 0; ci < (int) children.size(); ++ci)
	    branch_string[n] << '\n' << branch_string[children[ci]];
	  branch_string[n] << '\n' << nested_indent_string;
	}

      branch_string[n] << ')';
    }
  out << '\n' << branch_string[0] << '\n';

  // print constrained & unconstrained cliques
  if (path.size())
    {
      if (free_clique >= 0)
	{
	  out << base_indent_string;
	  show_clique (free_clique, false, out, show_composite_paths);
	}
      for (int n_clique = 0; n_clique < (int) clique.size(); ++n_clique)
	if (n_clique != free_clique)
	  {
	    out << base_indent_string;
	    show_clique (n_clique, true, out, show_composite_paths);
	  }
    }

  // print Redelings-Suchard proposal schedule
  if (proposal_branches.size())
    {
      out << base_indent_string << '(' << TSEXPR_SCHEDULE;
      for_const_contents (RedSuchSchedule, proposal_branches, pb)
	{
	  out << '\n' << base_indent_string << " (";
	  for (int i = 0; i < (int) pb->size(); ++i)
	    out << (i == 0 ? "" : " ") << branch_name[(*pb)[i]];
	  out << ')';
	}
      out << ")\n";
    }
}

void Transducer_SExpr_file::show_path_with_tree (const EHMM_transducer_funcs& ehmm, const vector<int>& path, Loge path_ll, ostream& out, int base_indent)
{
  // build tree string
  vector<sstring> branch_string (etree.nodes());
  const vector<vector<int> > branch_path = ehmm.branch_paths (path);
  NodePathMap branch_path_map;
  for (int n = 0; n < (int) branch_path.size(); ++n)
    branch_path_map[n] = branch_path[n];

  // output
  out << sstring (base_indent, ' ');

  show_path (ehmm, path, TSEXPR_STATE, out);
  out << '\n' << sstring (base_indent, ' ');

  show_path_types (ehmm, path, out);
  out << '\n';

  show_path_with_tree (branch_path_map, path_ll, out, base_indent);
}

void Transducer_SExpr_file::show_path_with_tree (const NodePathMap& branch_path, Loge path_ll, ostream& out, int base_indent)
{
  // build tree string
  vector<sstring> branch_string (etree.nodes());

  if (branch_path.size())
    for (int n = etree.nodes() - 1; n >= 0; --n)
      {
	int nested_indent = base_indent + 1;  // the extra +1 is so that the tree breakdown itself is indented
	for (int a = n; a > 0; a = etree.parent[a])
	  ++nested_indent;
	const sstring nested_indent_string (nested_indent, ' ');

	branch_string[n] << '\n' << nested_indent_string
			 << '(' << TSEXPR_BRANCH
			 << " (" << TSEXPR_NAME << ' ' << branch_name[n]
			 << ")\n" << nested_indent_string << ' ';

	show_path (branch_trans[n], ((NodePathMap&)branch_path)[n], TSEXPR_STATE, branch_string[n]);
	branch_string[n] << '\n' << nested_indent_string << ' ';

	show_path_types (branch_trans[n], ((NodePathMap&)branch_path)[n], branch_string[n]);

	const vector<ENode> children = etree.children (n);
	if (children.size())
	  branch_string[n] << '\n' << nested_indent_string << ' ';

	for (int ci = 0; ci < (int) children.size(); ++ci)
	  branch_string[n] << branch_string[children[ci]];

	branch_string[n] << ')';
      }

  // output
  out << sstring (base_indent, ' ');
  out << '(' << TSEXPR_FULL_SCORE << ' ' << -Nats2Bits (path_ll) << ')';

  out << '\n' << sstring (base_indent, ' ');
  out << branch_string[0];
}

void Transducer_SExpr_file::show_composite_name_format (ostream& out)
{
  vector<sstring> format_string (etree.nodes());
  for (int n = etree.nodes() - 1; n >= 0; --n)
    {
      format_string[n] << branch_name[n];

      const vector<ENode> children = etree.children (n);
      if (children.size())
	{
	  for (int ci = 0; ci < (int) children.size(); ++ci)
	    format_string[n] << (ci > 0 ? ' ' : '(') << format_string[children[ci]];
	  format_string[n] << ')';
	}
    }

  out << '(' << TSEXPR_FORMAT << " (" << format_string[0] << "))\n";
}

Score_profile Transducer_SExpr_file::peel_clique (int n_clique, const vector<Pair_transducer_scores>& branch_sc, Loge& clique_loglike, int prof_node, bool normalize_prof_node)
{
  // get clique info
  set<int>& clique_set = clique[n_clique];

  vector<bool> in_clique (etree.nodes(), false);
  vector<int> clique_vec;
  vector<int> path_cursor (etree.nodes(), 0);
  vector<int> seq_cursor (etree.nodes(), 0);
  for_const_contents (set<int>, clique_set, n)
    {
      in_clique[*n] = true;
      clique_vec.push_back (*n);
      if (prof_node < 0 && hinge_nodes.find (*n) != hinge_nodes.end())
	prof_node = *n;
    }
  const int clique_root = clique_vec[0];
  const bool root_is_hinged = hinge_nodes.find (clique_root) != hinge_nodes.end();

  // check that clique has something to peel
  if (clique_set.size() < 1 || (clique_set.size() == 1 && clique_root != 0))
    THROWEXPR ("Attempt to peel a clique that has no branches");

  // make temporary copy of branch_sc and (if clique_root is "hinged") place a singleton transducer at clique_root
  // (actually, since singleton is the most lightweight transducer, we put singletons at every node, then remove them selectively)
  Singleton_transducer_scores singleton (alphabet_size);
  vector<Pair_transducer_scores> tmp_branch_sc (branch_sc.size(), singleton);
  for_const_contents (vector<int>, clique_vec, n)
    if (*n != clique_root || !root_is_hinged)
      tmp_branch_sc[*n] = branch_sc[*n];

  // create peeler
  Transducer_peeler peeler;
  peeler.init (alphabet_size, etree, tmp_branch_sc);

  // infer length of clique root sequence
  int clique_root_length;
  if (prof.find (clique_root) != prof.end())
    clique_root_length = prof[clique_root].size();
  else
    {
      const vector<ENode> root_kids = etree.children (clique_root);
      int first_clique_root_child = -1;
      for_const_contents (vector<ENode>, root_kids, c)
	if (clique_set.find (*c) != clique_set.end())
	  {
	    first_clique_root_child = *c;
	    break;
	  }

      if (first_clique_root_child >= 0)
	{
	  clique_root_length = 0;
	  for_const_contents (vector<int>, path[first_clique_root_child], s)
	    {
	      const int t = branch_sc[first_clique_root_child].state_type[*s];
	      if (t == TransducerMatchType || t == TransducerDeleteType)
		++clique_root_length;
	    }
	}
      else
	THROWEXPR ("Weird -- clique root lacks any children -- this shouldn't be happening, man!");
    }
  
  // make temporary copy of path; if clique_root is hinged, place a dummy path of clique_root_length insert states at clique_root
  NodePathMap tmp_path;   // this is a bit inefficient
  for_const_contents (vector<int>, clique_vec, n)
    if (*n != clique_root || !root_is_hinged)
      tmp_path[*n] = path[*n];

  if (root_is_hinged)
    tmp_path[clique_root] = vector<int> (clique_root_length, singleton.insert_state());

  // debug: output tmp_path
  if (CTAGGING(2,TRANSDUCER_PEEL_DEBUG))
    for_const_contents (NodePathMap, tmp_path, np)
      CL << "tmp_path[" << np->first << "] = (" << np->second << ")\n";

  // align branch paths, do peeling, get posterior distribution at hinge node, calculate clique score
  vector<int> branch_state (etree.nodes(), (int) Start);
  vector<const Symbol_score_map*> node_scores (etree.nodes(), (const Symbol_score_map*) 0);

  Composite_path comp_path;
  comp_path.tmp_branch_sc = tmp_branch_sc;

  Score_profile result;
  clique_loglike = (Loge) 0;
  while (1)
    {
      // debug
      if (CTAGGING(0,TRANSDUCER_PEEL_DEBUG))
	CL << "branch_state=(" << branch_state << ") path_cursor=(" << path_cursor << ")\n";

      // find next branch transducer to advance: highest numbered insert state
      ENode inserter = -1;
      bool found_ins = false;
      for_const_reverse_contents (vector<int>, clique_vec, n)
	{
	  if (tmp_path.find(*n) == tmp_path.end())
	    THROWEXPR ("Attempt to advance cursor for nonexistent state path (tape " << tape_name[*n] << ")");
	  if (path_cursor[*n] < (int) tmp_path[*n].size())
	    if (tmp_branch_sc[*n].state_type[tmp_path[*n][path_cursor[*n]]] == Transducer_state_type_enum::TransducerInsertType)
	      {
		inserter = *n;
		found_ins = true;
		break;
	      }
	}
      if (!found_ins)
	break;

      // advance branches, propagating emissions downwards within the clique
      vector<ENode> deleters, matchers;
      stack<int> branches_to_advance;
      branches_to_advance.push (inserter);
      bool emission_to_prof_node = false;  // will be set if this column is relevant to the node we're trying to profile
      Loge trans_loglike = 0.;
      Composite_path_step comp_path_step;
      while (!branches_to_advance.empty())
	{
	  const int n = branches_to_advance.top();
	  if (path_cursor[n] >= (int) tmp_path[n].size())
	    {
	      show_tree (CLOGERR);
	      THROWEXPR ("Specified state path for branch " << branch_name[n] << " is inconsistent (too short?)");
	    }

	  const int s = tmp_path[n][path_cursor[n]];
	  const int t = tmp_branch_sc[n].state_type[s];

	  const Loge branch_trans_loglike = Score2Nats (tmp_branch_sc[n].effective_trans_score (branch_state[n], s));
	  NatsPMulAcc (trans_loglike, branch_trans_loglike);
	  if (CTAGGING(2,TRANSDUCER_PEEL_DEBUG))
	    CL << "Adding " << Nats2Bits(branch_trans_loglike) << " bits to transition log-likelihood (branch " << tape_name[n] << ")\n";

	  branch_state[n] = s;
	  branches_to_advance.pop();
	  ++path_cursor[n];

	  // update comp_path_step
	  comp_path_step.branch_trans[n] = s;

	  // is there an emission at this node?
	  if (t == Transducer_state_type_enum::TransducerInsertType || t == Transducer_state_type_enum::TransducerMatchType)
	    {
	      // if this is the node we're trying to profile, then set the appropriate flag
	      if (n == prof_node)
		emission_to_prof_node = true;

	      // advance node cursor and populate node_scores, if appropriate
	      if (prof.find (n) != prof.end())
		{
		  if (seq_cursor[n] >= (int) prof[n].size())
		    THROWEXPR ("Specified sequence for tape " << tape_name[n] << " is inconsistent (too short?)");

		  node_scores[n] = &prof[n][seq_cursor[n]];
		  ++seq_cursor[n];
		}

	      // put kids on "advance branch" stack
	      const vector<ENode> kids = etree.children (n);
	      for_const_reverse_contents (vector<ENode>, kids, c)
		if (clique_set.find (*c) != clique_set.end())
		  branches_to_advance.push (*c);
	    }

	  // init matchers, deleters
	  if (t == Transducer_state_type_enum::TransducerMatchType)
	    matchers.push_back (n);
	  else if (t == Transducer_state_type_enum::TransducerDeleteType)
	    deleters.push_back (n);

	  // update comp_path_step
	  comp_path_step.emit_nodes.push_back (n);
	}

      // ensure matchers & deleters are in the right order
      reverse (matchers.begin(), matchers.end());
      reverse (deleters.begin(), deleters.end());

      // log
      CTAG(2,TRANSDUCER TRANSDUCER_PEEL) << "At position " << result.size() << " of tape " << prof_node << " (" << tape_name[prof_node]
					 << "): path_cursor = (" << path_cursor << "), seq_cursor = (" << seq_cursor << ")\n";

      // create emission; prune
      Transducer_tree_emission ttemit (node_scores, branch_state, inserter, deleters, matchers);
      peeler.prune (ttemit);

      // count transition & emission in cumulative log-likelihood
      const bool ignore_emission = emission_to_prof_node && !normalize_prof_node;
      NatsPMulAcc (clique_loglike, NatsPMul (trans_loglike, (ignore_emission ? 0. : Score2Nats (peeler.pruning_sc))));
      CTAG(2,TRANSDUCER TRANSDUCER_PEEL) << "Transition " << Nats2Bits(trans_loglike) << " bits; emission " << Score2Bits (peeler.pruning_sc) << " bits"
					 << (ignore_emission ? " (emission ignored)" : "") << '\n';

      // emission to profile node?
      if (emission_to_prof_node)
	{
	  // peel
	  peeler.peel (ttemit);

	  // get posterior distribution; push onto result profile
	  result.push_back (peeler.post_sc (prof_node, normalize_prof_node));
	  if (CTAGGING(2,TRANSDUCER TRANSDUCER_PEEL))
	    {
	      CL << "Profile: ";
	      show_ssm (result.back(), CL);
	      CL << '\n';
	    }
	}
      else
	CTAG(2,TRANSDUCER TRANSDUCER_PEEL) << "(nothing to see here; moving along)\n";

      // update comp_path
      comp_path.steps.push_back (comp_path_step);
    }

  // check that all cursors reached the end of their respective paths, and count end transitions in cumulative score
  Loge final_trans_loglike = 0.;
  for_const_contents (vector<int>, clique_vec, n)
    {
      if (prof.find (*n) != prof.end())
	if (seq_cursor[*n] < (int) prof[*n].size())
	  {
	    show_tree (CLOGERR);
	    THROWEXPR ("Specified sequence for tape " << tape_name[*n] << " is inconsistent (too long?)");
	  }
      
      if (path_cursor[*n] < (int) tmp_path[*n].size())
	{
	  show_tree (CLOGERR);
	  THROWEXPR ("Specified state path for branch " << branch_name[*n] << " is inconsistent (too long?)");
	}

      NatsPMulAcc (final_trans_loglike, Score2Nats (tmp_branch_sc[*n].effective_trans_score (branch_state[*n], End)));
    }

  // store composite path for this clique
  clique_path[n_clique] = comp_path;

  // log final transition score
  CTAG(2,TRANSDUCER TRANSDUCER_PEEL) << "Final transition " << Nats2Bits(final_trans_loglike) << " bits\n";
  NatsPMulAcc (clique_loglike, final_trans_loglike);

  // log clique score
  CTAG(2,TRANSDUCER TRANSDUCER_PEEL) << "Clique log-likelihood is " << Nats2Bits(clique_loglike) << " bits\n";

  // return
  return result;
}

vector<Pair_transducer_scores> Transducer_SExpr_file::make_branch_sc()
{
  vector<Pair_transducer_scores> branch_sc;
  for_const_contents (vector<Pair_trans_funcs>, branch_trans, bt)
    branch_sc.push_back (bt->eval_sc (pscores));
  return branch_sc;
}

Transducer_SExpr_file Transducer_SExpr_file::peel_constrained (Loge& peeling_loglike, bool normalize_prof_nodes)
{
  if (free_clique < 0)
    THROWEXPR ("There is no subtree on whose branches state paths are not specified");

  const vector<Pair_transducer_scores> branch_sc = make_branch_sc();

  const vector<int> subtree_nodes (clique[free_clique].begin(), clique[free_clique].end());
  ETree subtree = etree.subtree (subtree_nodes);

  vector<Pair_trans_funcs> subtree_branch_trans;
  map<int,sstring> subtree_tape_name;
  vector<sstring> subtree_branch_trans_name, subtree_branch_name;
  NodeProfileMap subtree_prof;
  NodePathMap subtree_path, subtree_old_path;
  map<int,double> subtree_band_coeff;
  vector<int> tree2subtree (etree.nodes(), -1);

  for (int f = 0; f < (int) subtree_nodes.size(); ++f)
    {
      const int n = subtree_nodes[f];
      tree2subtree[n] = f;

      subtree_tape_name[f] = tape_name[n];
      subtree_branch_trans.push_back (branch_trans[n]);
      subtree_branch_name.push_back (branch_name[n]);
      subtree_branch_trans_name.push_back (branch_trans_name[n]);

      if (old_path.find (n) != old_path.end())
	subtree_old_path[f] = old_path[n];

      if (band_coeff.find (n) != band_coeff.end())
	subtree_band_coeff[f] = band_coeff[n];

      if (hinge_nodes.find (n) != hinge_nodes.end())
	{
	  Loge clique_ll = -InfinityLoge;
	  subtree_prof[f] = peel_clique (cons_clique[n], branch_sc, clique_ll, -1, normalize_prof_nodes);
	  NatsPMulAcc (peeling_loglike, clique_ll);
	}
      else if (prof.find (n) != prof.end())
	subtree_prof[f] = prof[n];
    }

  // if the branch leading to the free subtree's root is constrained (so that the original singleton transducer has been factored in during peeling),
  // place a dummy singleton transducer on the clique root branch; also add an eponymous singleton transducer to our dictionary, so that it's comprehensive.
  // NB root_is_hinged is NOT just the same as (clique_root != 0), since node 0 can be hinged too.
  const int clique_root = subtree_nodes[0];
  const bool root_branch_constrained = path.find (clique_root) != path.end();
  if (root_branch_constrained)
    {
      sstring singleton_name (Singleton_transducer_name);
      subtree_branch_trans_name[0] = singleton_name;
      Singleton_transducer_funcs sing_funcs (alphabet_size);
      trans_dict[singleton_name] = subtree_branch_trans[0] = sing_funcs;
      subtree_path[0] = vector<int> (subtree_prof[0].size(), sing_funcs.insert_state());
    }

  // to re-use root node name, comment out the next line of code and add the following:
  //  if (subtree_nodes.size())
  //    subtree_tape_name[-1] = tape_name[etree.parent[clique_root]];
  if (subtree_branch_name.size())
    subtree_branch_name[0].clear();

  // set up Redelings-Suchard schedule
  RedSuchSchedule subtree_proposal_branches;
  for_const_contents (RedSuchSchedule, proposal_branches, pb)
    {
      vector<int> subtree_branch;
      for_const_contents (vector<int>, *pb, b)
	subtree_branch.push_back (tree2subtree[*b]);
      subtree_proposal_branches.push_back (subtree_branch);
    }

  // build peeled Transducer_SExpr_file object
  Transducer_SExpr_file subtree_transducer_sexpr_file (alphabet, pscores, subtree, subtree_branch_trans, subtree_prof, subtree_path, subtree_old_path);
  subtree_transducer_sexpr_file.proposal_branches = subtree_proposal_branches;
  subtree_transducer_sexpr_file.band_coeff = subtree_band_coeff;

  // rebuild tree names
  subtree_transducer_sexpr_file.rebuild_tree_names (subtree_tape_name, subtree_branch_name, subtree_branch_trans_name);

  // return
  return subtree_transducer_sexpr_file;
}

void Transducer_SExpr_file::rebuild_tree_names (const map<int,sstring>& new_tape_name, const vector<sstring>& new_branch_name, const vector<sstring>& new_branch_trans_name)
{
  tape_name = new_tape_name;
  branch_name = new_branch_name;
  branch_trans_name = new_branch_trans_name;

  autoname_tree();
  autobuild_trans_dict();
  setup_composite_name();
}

void Transducer_SExpr_file::setup_cliques()
{
  free_clique = -1;
  uncons_clique.clear();
  cons_clique.clear();
  hinge_nodes.clear();
  branch_cons.clear();
  clique.clear();

  for (int node = 0; node < etree.nodes(); ++node)
    {
      const int parent = etree.parent[node];
      const bool node_has_path = path.find (node) != path.end();
      if (CTAGGING(-1,TRANSDUCER_DEBUG))
	CL << "Node " << node << " (parent " << parent << ") " << (node_has_path ? "has" : "does not have") << " a path\n";
      map<int,int>& node_clique = node_has_path ? cons_clique : uncons_clique;
      int n_clique = clique.size();
      if (parent >= 0)
	{
	  if (node_clique.find (parent) != node_clique.end())
	    n_clique = node_clique[parent];
	  else
	    {
	      node_clique[parent] = n_clique;
	      hinge_nodes.insert (parent);
	    }
	}
      node_clique[node] = n_clique;
      branch_cons.push_back (node_has_path);

      if (n_clique == (int) clique.size())
	clique.push_back (set<int>());
      clique[n_clique].insert (node);
      if (parent >= 0)
	clique[n_clique].insert (parent);

      if (!node_has_path)
	{
	  if (free_clique >= 0 && free_clique != n_clique)
	    THROWEXPR ("Can't have more than one disjoint subtree with no specified state path");
	  free_clique = n_clique;
	}
    }
}

void Transducer_SExpr_file::show_clique (int n_clique, bool constrained, ostream& out, bool show_composite_paths)
{
  set<int>& clique_set = clique[n_clique];
  const int clique_root = *clique_set.begin();

  vector<bool> in_clique (etree.nodes(), false);
  vector<int> clique_vec;
  sstring profile_string;
  for_const_contents (set<int>, clique_set, n)
    {
      in_clique[*n] = true;
      clique_vec.push_back (*n);
      if (constrained && hinge_nodes.find (*n) != hinge_nodes.end())
	profile_string << tape_name[*n];
    }

  vector<sstring> clique_string (etree.nodes());
  for_const_reverse_contents (vector<int>, clique_vec, n)
    {
      vector<int> kids_in_clique;
      const vector<int> kids = etree.children (*n);
      for_const_contents (vector<int>, kids, c)
	if (in_clique[*c])
	  kids_in_clique.push_back (*c);

      const bool hide_branch = *n == clique_root && branch_cons[*n] != constrained;
      if (!hide_branch)
	clique_string[*n] << branch_name[*n];

      if (kids_in_clique.size())
	{
	  if (!hide_branch)
	    clique_string[*n] << '(';
	  for (int i = 0; i < (int) kids_in_clique.size(); ++i)
	    {
	      if (i > 0) clique_string[*n] << ' ';
	      clique_string[*n] << clique_string[kids_in_clique[i]];
	    }
	  if (!hide_branch)
	    clique_string[*n] << ')';
	}
    }

  out << '(' << (constrained ? TSEXPR_CONS : TSEXPR_UNCONS);
  if (profile_string.size())
    out << "\n (" << TSEXPR_PROFILE << ' ' << profile_string << ')';
  out << "\n (" << TSEXPR_BRANCH << " (" << clique_string[clique_root] << "))";
  if (constrained && show_composite_paths && clique_path.find (n_clique) != clique_path.end())
    {
      out << "\n (" << TSEXPR_ALIGNMENT;
      const sstring indent_string (2, ' '); 
      Composite_path& comp_path = clique_path[n_clique];
      for (int i = 0; i < (int) comp_path.steps.size(); ++i)
	{
	  Composite_path_step& comp_path_step = comp_path.steps[i];
	  out << '\n' << indent_string
	      << '(' << TSEXPR_ALIGN_COL
	      << " (" << TSEXPR_ID << ' ' << i+1 << ')';
	  for_const_contents (Branch_transition_map, comp_path_step.branch_trans, bt)
	    {
	      out << "\n " << indent_string;
	      out << '(' << TSEXPR_TRANSITION
		  << " (" << TSEXPR_BRANCH << ' ' << branch_name[bt->first]
		  << ") (" << TSEXPR_TO << ' ' << comp_path.tmp_branch_sc[bt->first].get_state_name (bt->second)
		  << "))";
	    }
	  out << '\n' << indent_string << " (" << TSEXPR_TYPE << " (";

	  int n_emit_node = 0;
	  for_const_contents (vector<int>, comp_path_step.emit_nodes, n)
	    out << (n_emit_node++ ? " " : "") << tape_name[*n];
	  out << ")))";
	}
      out << ")";
   } 
  out << ")\n";
}


void Transducer_SExpr_file::show_composite (EHMM_transducer_funcs& ehmm, ostream& out)
{
  // do we have any PVar definitions?
  const bool show_bitvalues = defined_pvars.size() > 0;

  // print composite transducer
  out << '\n';
  ehmm.show_sexpr (out, composite_name, &pscores, show_bitvalues);
}

void Transducer_SExpr_file::show_prof (const Score_profile& prof_sc, ostream& out, int base_indent)
{
  bool is_bitmatrix = false;
  for_const_contents (Score_profile, prof_sc, ssm)
    if (ssm->size() == 1 ? ssm->begin()->second != 0 : true)
      {
	is_bitmatrix = true;
	break;
      }

  const char* tag = is_bitmatrix ? TSEXPR_BITPROFILE : TSEXPR_SEQUENCE;

  out << '(' << tag << " (";;

  sstring spacer;
  if (is_bitmatrix)
    spacer << '\n' << sstring (base_indent + strlen(tag) + 3, ' ');
  else
    spacer << ' ';

  int pos = 0;
  for_const_contents (Score_profile, prof_sc, ssm)
    {
      if (pos++)
	out << spacer;
      if (ssm->size() == 1 && ssm->begin()->second == 0)
	out << alphabet[ssm->begin()->first];
      else
	{
	  out << '(';
	  int n_sym = 0;
	  for_const_contents (Symbol_score_map, *ssm, ss)
	    if (ss->second > -InfinityScore)
	      {
		if (n_sym++)
		  out << ' ';
		out << '(' << alphabet[ss->first] << ' ' << score_sexpr(ss->second) << ')';
	      }
	  out << ')';
	}
    }
  out << "))";
}

void Transducer_SExpr_file::show_ssm (const Symbol_score_map& ssm, ostream& out)
{
  if (ssm.size() == 1 && ssm.begin()->second == 0)
    out << alphabet[ssm.begin()->first];
  else
    {
      out << '(';
      for_const_contents (Symbol_score_map, ssm, ss)
	if (ss->second > -InfinityScore)
	  out << '(' << alphabet[ss->first] << ' ' << score_sexpr(ss->second) << ')';
      out << ')';
    }
}

vector<Score_profile*> Transducer_SExpr_file::make_seq_vec()
{
  vector<Score_profile*> seq (etree.nodes(), (Score_profile*) 0);
  for_contents (NodeProfileMap, prof, node_prof)
    {
      seq[node_prof->first] = &node_prof->second;
      if (CTAGGING(0,TRANSDUCER TRANSDUCER_SEQVEC))
	{
	  CL << "make_seq_vec[" << node_prof->first << "]:\n";
	  show_prof (*seq[node_prof->first], CL, 0);
	  CL << '\n';
	}
    }
  return seq;
}

vector<int> Transducer_SExpr_file::observed_nodes()
{
  vector<int> obs;
  for_contents (NodeProfileMap, prof, node_prof)
    obs.push_back (node_prof->first);
  return obs;
}

sstring Transducer_SExpr_file::score_sexpr (Score sc)
{
  return PFunc_builder::score_sexpr (sc, true);
}

sstring Transducer_SExpr_file::score_sexpr (Loge ll)
{
  sstring s;
  if (ll <= -InfinityLoge)
    s << PK_INFINITE;
  else
    {
      s.precision(3);
      s << fixed << -Nats2Bits (ll);
    }
  return s;
}

Singleton_transducer_funcs::Singleton_transducer_funcs (int alph_sz)
  : Pair_transducer_funcs (2)
{
  Pair_transducer_funcs::alphabet_size = alph_sz;

  state_type[insert_state()] = TransducerInsertType;
  state_type[wait_state()] = TransducerWaitType;
  alloc_pair_emit (PFunc (1));

  start_name = Singleton_start_state_name;
  end_name = Singleton_end_state_name;
  state_name[insert_state()] = Singleton_insert_state_name;
  state_name[wait_state()] = Singleton_wait_state_name;

  transition (Start, insert_state()) = 1;
  transition (Start, wait_state()) = 1;

  transition (insert_state(), insert_state()) = 1;
  transition (insert_state(), wait_state()) = 1;

  transition (wait_state(), End) = 1;
}

Singleton_transducer_scores::Singleton_transducer_scores (int alph_sz)
  : Pair_transducer_scores (2)
{
  Pair_transducer_scores::alphabet_size = alph_sz;

  state_type[insert_state()] = TransducerInsertType;
  state_type[wait_state()] = TransducerWaitType;
  alloc_pair_emit (0);

  start_name = Singleton_start_state_name;
  end_name = Singleton_end_state_name;
  state_name[insert_state()] = Singleton_insert_state_name;
  state_name[wait_state()] = Singleton_wait_state_name;

  transition (Start, insert_state()) = 0;
  transition (Start, wait_state()) = 0;

  transition (insert_state(), insert_state()) = 0;
  transition (insert_state(), wait_state()) = 0;

  transition (wait_state(), End) = 0;
}

void Transducer_SExpr_file::simulate()
{
  if (prof.size() || path.size())
    THROWEXPR ("Can't sample paths/sequences, because some were already specified");
  Digitized_biosequence empty_dsq;
  map<int,Digitized_biosequence> dsq;
  for (ENode node = 0; node < etree.nodes(); ++node)
    {
      const Digitized_biosequence& parent_dsq (node == 0 ? empty_dsq : dsq[etree.parent[node]]);
      vector<int> parent_child_path;
      Digitized_biosequence child_dsq;
      const Pair_transducer_funcs& parent_child_trans_funcs (branch_trans[node]);
      Pair_transducer_scores parent_child_trans_sc = parent_child_trans_funcs.eval_sc (pscores);
      parent_child_trans_sc.sample (parent_dsq, parent_child_path, child_dsq);
      dsq[node] = child_dsq;
      prof[node] = Score_profile (child_dsq);
      path[node] = parent_child_path;
    }
}


Stockade Transducer_SExpr_file::stockade()
{
  // make the alphabet "flush", so all symbols are the same length
  int sym_len = 0;
  for_const_contents (vector<sstring>, alphabet, sym)
    if ((int) sym->size() > sym_len)
      sym_len = sym->size();

  vector<sstring> flush_alph;
  for_const_contents (vector<sstring>, alphabet, sym)
    {
      sstring flush_sym (*sym);
      while ((int) flush_sym.size() < sym_len)
	flush_sym << Alignment::gap_char();
      flush_alph.push_back (flush_sym);
    }

  // create Alignment_path (via Alignment_path::Decomposition)
  Alignment_path::Decomposition decomp;
  for (int node = 1; node < etree.nodes(); ++node)
    {
      Pairwise_path pair_path;
      if (path.find (node) == path.end())
	{
	  // no path found, so assume an all-gap alignment
	  if (prof.find (etree.parent[node]) != prof.end())
	    {
	      const int parent_len = prof[etree.parent[node]].size();
	      for (int n = 0; n < parent_len; ++n)
		pair_path.append_column (1, 0);
	    }

	  if (prof.find (node) != prof.end())
	    {
	      const int child_len = prof[node].size();
	      for (int n = 0; n < child_len; ++n)
		pair_path.append_column (0, 1);
	    }
	}
      else
	{
	  // path found, so convert state path into pairwise alignment
	  pair_path = branch_trans[node].state2align_path (path[node], true);
	}

      const Alignment_path::Row_pair row_pair (etree.parent[node], node);
      decomp[row_pair] = pair_path;
    }

  // build the multi-alignment path out of the decomposition
  Alignment_path align_path;
  align_path.compose_and_log (decomp, false, true);

  // create Stockade
  Stockade stockade (align_path.rows(), align_path.columns());
  stockade.align.path = align_path;

  // create PHYLIP_tree
  PHYLIP_tree tree;
  for (int node = 0; node < etree.nodes(); ++node)
    {
      tree.add_node (etree.parent[node]);
      tree.node_name.push_back (tape_name[node]);
    }
  tree.root = 0;
  tree.rebuild_parents();

  // add tree annotation to Stockholm alignment
  sstring tree_string;
  tree.write (tree_string, 0);
  stockade.align.add_gf_annot (sstring (Stockholm_New_Hampshire_tag), tree_string);

  // set row names
  for (int node = 0; node < etree.nodes(); ++node)
    stockade.align.row_name[node] = tape_name[node];

  // populate Biosequence's in Named_profiles
  for (int node = 0; node < etree.nodes(); ++node)
    {
      if (prof.find (node) != prof.end())
	{
	  const Score_profile& prof_sc = prof[node];
	  for (int pos = 0; pos < (int) prof_sc.size(); ++pos)
	    stockade.np[node].seq << flush_alph [prof_sc.consensus (pos)];
	}
      else
	{
	  // set np[node] and prof[node] to null pointers
	  // this will display Felsenstein wildcards instead
	  stockade.align.prof[node] = (Score_profile*) 0;
	  stockade.align.np[node] = (Named_profile*) 0;
	}
    }

  // return
  return stockade;
}

Loge Transducer_SExpr_file::get_path_loglike (const NodePathMap& path_map) const
{
  Loge ll;

  NodePathMap dummy_path_map, full_path_map (old_path);
  for_const_contents (NodePathMap, path_map, pm)
    full_path_map[pm->first] = pm->second;

  if ((int) full_path_map.size() < etree.nodes() - 1)
    THROWEXPR ("Path underspecified");

  Transducer_SExpr_file path_comp (alphabet, pscores, etree, branch_trans, prof, full_path_map, dummy_path_map);
  path_comp.rebuild_tree_names (tape_name, branch_name, branch_trans_name);

  const vector<Pair_transducer_scores> branch_sc = path_comp.make_branch_sc();
  path_comp.peel_clique (0, branch_sc, ll, 0);

  return ll;
}

Transducer_SExpr_file::NodePathMap Transducer_SExpr_file::get_path_map (const EHMM_transducer_funcs& ehmm, const vector<int>& trace) const
{
  const vector<vector<int> > branch_path = ehmm.branch_paths (trace);
  
  Transducer_SExpr_file::NodePathMap path_map;
  for (int n = 0; n < ehmm.etree.nodes(); ++n)
    path_map[n] = branch_path[n];

  return path_map;
}

void Transducer_SExpr_file::autopropose()
{
  // outline of algorithm:
  // - let H be the set of all hinge nodes, or nodes for which sequence is specified
  // - let S be H minus the clique root
  // - let K be the set of "known" nodes, initially equal to H-S
  // - while S is nonempty:
  //  - let N be [the first] node in S
  //  - if N has any ancestors in K:
  //   - let A be N's most recent ancestor in K
  //   - add path from A-->N to proposal schedule
  //  - else (N has no ancestors in K):
  //   - let M be [the second] node in S
  //   - let A be the MRCA of M,N
  //   - add paths from A-->M and A-->N to proposal schedule
  //   - if A is not the clique root, add A to S
  //   - if A==0 (i.e. A is the global root) and the subroot->root path is unconstrained, add subroot->root to proposal schedule
  //   - remove M from S
  //  - remove N from S
  //  - add all nodes from current proposal to K

  proposal_branches.clear();

  if (free_clique >= 0)
    {
      const set<int>& f = clique[free_clique];

      // initialize H
      set<int> H = hinge_nodes;
      for_const_contents (NodeProfileMap, prof, np)
	if (f.find (np->first) != f.end())
	  H.insert (np->first);

      // initialize S and K
      set<int> S (H), K;
      // remove clique root from S, if it's there
      int s_min = *(S.begin());  // if it's in S, clique root will be minimum-valued node in S
      if (f.find (etree.parent[s_min]) == f.end())  // the clique root's parent is not in S
	{
	  S.erase (s_min);
	  K.insert (s_min);
	}

      // main loop
      vector<int> branches;
      while (S.size())
	{
	  set<int> current_path;
	  set<int>::const_iterator S_iter = S.begin();
	  const int N = *S_iter;
	  int n;
	  for (n = N; n >= 0 && K.find (n) == K.end(); n = etree.parent[n])
	    current_path.insert (n);
	  if (n < 0)
	    {
	      set<int> N_ancestors;
	      swap (N_ancestors, current_path);
	      const int M = *(++S_iter);
	      int m;
	      for (m = M; N_ancestors.find (m) == N_ancestors.end(); m = etree.parent[m])
		current_path.insert (m);
	      for (n = N; n != m; n = etree.parent[n])
		current_path.insert (n);
	      // at this point m=A; and if m's parent is in f, then m is not the clique root
	      if (f.find (etree.parent[m]) != f.end())
		S.insert (m);
	      // if m==0 (i.e. it's the global root) and there's no path from the subroot, then add m to path
	      if (m == 0 && path.find (0) == path.end())
		current_path.insert (0);
	    }
	  for_const_contents (set<int>, current_path, b)
	    {
	      S.erase (*b);
	      K.insert (*b);
	    }
	  branches.insert (branches.end(), current_path.begin(), current_path.end());
	  proposal_branches.push_back (branches);
	}
    }
}

bool Transducer_SExpr_file::has_all_banding_coefficients() const
{
  const set<int>& fcs = free_clique_set();
  for_const_contents (set<int>, fcs, n)
    {
      const int p = etree.parent[*n];
      if (p >= 0 && fcs.find (p) != fcs.end())  // exclude root of clique; we're interested in branches
	if (band_coeff.find (*n) == band_coeff.end())
	  return false;
    }
  return true;
}

long Transducer_SExpr_file::get_band_diameter() const
{
  // formula is: band_diameter = effective_tree_banding_coefficient * sqrt(rms_sequence_length)
  // where effective_tree_banding_coefficient = max_{leaf nodes L} sum_{nodes N from root to L} banding_coefficient[N]
  long band_diameter = 0;
  const set<int>& fcs = free_clique_set();
  if (fcs.size())
    {
      // compute effective_tree_banding_coefficient
      vector<double> band_coeff_sum (etree.nodes(), 0.);
      for_const_reverse_contents (set<int>, fcs, n)
	{
	  double max_child_coeff = 0.;
	  const vector<ENode> kids = etree.children (*n);
	  for_const_contents (vector<ENode>, kids, c)
	    if (fcs.find (*c) != fcs.end())
	      max_child_coeff = max (max_child_coeff, band_coeff_sum[*c]);

	  band_coeff_sum[*n] = max_child_coeff;
	  map<int,double>::const_iterator band_coeff_iter = band_coeff.find (*n);
	  if (band_coeff_iter != band_coeff.end())
	    band_coeff_sum[*n] += band_coeff_iter->second;
	}

      // compute rms_sequence_length
      double n_seq = 0., len_sq_sum = 0.;
      for_const_contents (NodeProfileMap, prof, node_prof)
	{
	  const int prof_size = node_prof->second.size();
	  n_seq += 1.;
	  len_sq_sum += prof_size * prof_size;
	}
      const double rms_len = n_seq == 0. ? 0. : sqrt (len_sq_sum / n_seq);

      // compute band_diameter
      const int clique_root = *(fcs.begin());
      band_diameter = (long) (1. + band_coeff_sum[clique_root] * sqrt (rms_len));  // round up
    }

  // return
  return band_diameter;
}
