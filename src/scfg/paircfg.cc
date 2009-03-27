#include <stack>
#include "scfg/paircfg.h"
#include "seq/stockholm.h"
#include "util/vector_output.h"

#define _CFGIDX_N 0
#define _CFGIDX_0 1
#define _CFGIDX_1 CFG_alphabet_size
#define _CFGIDX_2 (CFG_alphabet_size * CFG_alphabet_size)
#define _CFGIDX_3 (CFG_alphabet_size * CFG_alphabet_size * CFG_alphabet_size)
#define _CFGIDX_4 (CFG_alphabet_size * CFG_alphabet_size * CFG_alphabet_size * CFG_alphabet_size)

int Pair_CFG_state_typing::_emit_sz[EmitStateTypes] =
{ _CFGIDX_N, _CFGIDX_1, _CFGIDX_1, _CFGIDX_2,
  _CFGIDX_1, _CFGIDX_2, _CFGIDX_2, _CFGIDX_3,
  _CFGIDX_1, _CFGIDX_2, _CFGIDX_2, _CFGIDX_3,
  _CFGIDX_2, _CFGIDX_3, _CFGIDX_3, _CFGIDX_4 };

int Pair_CFG_state_typing::_emit_xl_mul[EmitStateTypes] =
{ _CFGIDX_N, _CFGIDX_0, _CFGIDX_N, _CFGIDX_0,
  _CFGIDX_N, _CFGIDX_0, _CFGIDX_N, _CFGIDX_0,
  _CFGIDX_N, _CFGIDX_0, _CFGIDX_N, _CFGIDX_0,
  _CFGIDX_N, _CFGIDX_0, _CFGIDX_N, _CFGIDX_0 };

int Pair_CFG_state_typing::_emit_xr_mul[EmitStateTypes] =
{ _CFGIDX_N, _CFGIDX_N, _CFGIDX_0, _CFGIDX_1,
  _CFGIDX_N, _CFGIDX_N, _CFGIDX_0, _CFGIDX_1,
  _CFGIDX_N, _CFGIDX_N, _CFGIDX_0, _CFGIDX_1,
  _CFGIDX_N, _CFGIDX_N, _CFGIDX_0, _CFGIDX_1 };

int Pair_CFG_state_typing::_emit_yl_mul[EmitStateTypes] =
{ _CFGIDX_N, _CFGIDX_N, _CFGIDX_N, _CFGIDX_N,
  _CFGIDX_0, _CFGIDX_1, _CFGIDX_1, _CFGIDX_2,
  _CFGIDX_N, _CFGIDX_N, _CFGIDX_N, _CFGIDX_N,
  _CFGIDX_0, _CFGIDX_1, _CFGIDX_1, _CFGIDX_2 };

int Pair_CFG_state_typing::_emit_yr_mul[EmitStateTypes] =
{ _CFGIDX_N, _CFGIDX_N, _CFGIDX_N, _CFGIDX_N,
  _CFGIDX_N, _CFGIDX_N, _CFGIDX_N, _CFGIDX_N,
  _CFGIDX_0, _CFGIDX_1, _CFGIDX_1, _CFGIDX_2,
  _CFGIDX_1, _CFGIDX_2, _CFGIDX_2, _CFGIDX_3 };

const char* Pair_CFG_state_type_enum::state_type_string (State_type t)
{
  switch (t)
    {
    case Null: return "Null"; break;
    case EmitXL: return "EmitXL"; break;
    case EmitXR: return "EmitXR"; break;
    case EmitXLR: return "EmitXLR"; break;
    case EmitYL: return "EmitYL"; break;
    case EmitYR: return "EmitYR"; break;
    case EmitYLR: return "EmitYLR"; break;
    case EmitXLYL: return "EmitXLYL"; break;
    case EmitXLYR: return "EmitXLYR"; break;
    case EmitXLYLR: return "EmitXLYLR"; break;
    case EmitXRYL: return "EmitXRYL"; break;
    case EmitXRYR: return "EmitXRYR"; break;
    case EmitXRYLR: return "EmitXRYLR"; break;
    case EmitXLRYL: return "EmitXLRYL"; break;
    case EmitXLRYR: return "EmitXLRYR"; break;
    case EmitXLRYLR: return "EmitXLRYLR"; break;
    case Bifurc: return "Bifurc"; break;
    case BifurcRevY: return "BifurcRevY"; break;
    case Undefined: return "Undefined"; break;
    default: return "[unknown]"; break;
    }
  return "UnreachableText";
}

void Pair_CFG_state_typing::init_bifurc (int state, int l, int r)
{
  state_type[state] = Bifurc;
  bifurc[state] = Bifurcation (l, r);
}

Pair_CFG_branch::Pair_CFG_branch (int xl, int xr, int yl, int yr, int p)
  : xl(xl), xr(xr), yl(yl), yr(yr), parent(p), lchild(-1), rchild(-1), revy(0)
{ }

int Pair_CFG_parse_tree::new_branch (int xl, int xr, int yl, int yr, int p)
{
  push_back (Pair_CFG_branch (xl, xr, yl, yr, p));
  return size() - 1;
}

int Pair_CFG_parse_tree::new_lchild (int xl, int xr, int yl, int yr, int p)
{
  const int lc = new_branch (xl, xr, yl, yr, p);
  (*this)[p].lchild = lc;
  return lc;
}

int Pair_CFG_parse_tree::new_rchild (int xl, int xr, int yl, int yr, int p)
{
  const int rc = new_branch (xl, xr, yl, yr, p);
  (*this)[p].rchild = rc;
  return rc;
}

bool Pair_CFG_parse_tree::test_connections() const
{
  if (size())
    {
      const Pair_CFG_branch& root = front();
      if (root.parent != -1) return 0;
      int n_children = 0;
      for (int i = 0; i < (int) size(); ++i)
	{
	  const Pair_CFG_branch& branch = (*this)[i];
	  if (branch.xl > branch.xr || branch.yl > branch.yr) return 0;
	  if (branch.path.size() == 0) return 0;  // every branch must have at least one state
	  if (branch.lchild >= 0)
	    {
	      if (branch.rchild < 0) return 0;
	      if (branch.lchild >= (int) size() || branch.rchild >= (int) size()) return 0;
	      if (branch.lchild == branch.rchild) return 0;
	      const Pair_CFG_branch& lchild = (*this)[branch.lchild];
	      const Pair_CFG_branch& rchild = (*this)[branch.rchild];
	      if (lchild.parent != i || rchild.parent != i) return 0;
	      if (lchild.xr != rchild.xl) return 0;
	      if (branch.revy ? (lchild.yl != rchild.yr) : (lchild.yr != rchild.yl)) return 0;
	      n_children += 2;
	    }
	  else
	    {
	      if (branch.rchild >= 0) return 0;
	      if (branch.revy) return 0;  // revy shouldn't be set on terminating branches
	    }
	}
      if (n_children != (int) size() - 1) return 0;
    }
  return 1;
}

bool Pair_CFG_parse_tree::test_global (const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy) const
{
  if (size() == 0) return 0;
  const Pair_CFG_branch& root = front();
  if (root.path.size() == 0) return 0;
  if (root.path[0] != HMM_state_enum::Start) return 0;
  if (root.xl != 0 || root.yl != 0 || root.xr != (int) dsqx.size() || root.yr != (int) dsqy.size()) return 0;
  return 1;
}

void Pair_CFG_parse_tree::show (ostream& out, const vector<State_type>* state_type) const
{
  typedef pair<int,int> Branch_indent;
  stack<Branch_indent> next_branch;
  next_branch.push (Branch_indent (0, 0));
  while (!next_branch.empty())
    {
      const Branch_indent bi = next_branch.top();
      const Pair_CFG_branch& b = (*this) [bi.first];
      const int indent = bi.second;
      next_branch.pop();
      for (int i = 0; i < indent; ++i)
	out << ' ';
      out << "[" << b.xl << ".." << b.xr << "," << b.yl << ".." << b.yr << "]";
      for_const_contents (vector<int>, b.path, s)
	{
	  out << " " << *s;
	  if (*s == Grammar_state_enum::Start)
	    out << "(Start)";
	  else if (*s == Grammar_state_enum::End)
	    out << "(End)";
	  else if (state_type)
	    out << "(" << state_type_string ((*state_type)[*s]) << ")";
	}
      out << "\n";
      if (b.has_children())
	{
	  next_branch.push (Branch_indent (b.rchild, indent+1));
	  next_branch.push (Branch_indent (b.lchild, indent+1));
	}
    }
}

void Pair_CFG_local_path::show (ostream& out) const
{
  out << "xstart=" << xstart << " xlen=" << xlen << " ystart=" << ystart << " ylen=" << ylen << " path=(" << path << ")\n";
}

Pair_CFG_local_path::Pair_CFG_local_path (int xstart, int xlen, int ystart, int ylen, const vector<int>& path)
  : xstart (xstart), xlen (xlen), ystart (ystart), ylen (ylen), path (path)
{ }

Pair_CFG_local_path::Pair_CFG_local_path()
{ }

Pair_CFG_state_typing& Pair_CFG_state_typing::operator= (const Pair_CFG_state_typing& base)
{
  assign_state_typing (base);
  return *this;
}

void Pair_CFG_state_typing::assign_state_typing (const Pair_CFG_state_typing& base)
{
  state_type = base.state_type;
  bifurc = base.bifurc;
}

void Pair_CFG_state_typing::indent_xr (Pair_CFG_parse_tree& parse_tree, int node, int offset)
{
  if (node < 0) return;
  parse_tree[node].xr += offset;
  indent_xr (parse_tree, parse_tree[node].lchild, offset);
  indent_xr (parse_tree, parse_tree[node].rchild, offset);
}

void Pair_CFG_state_typing::indent_yr (Pair_CFG_parse_tree& parse_tree, int node, int offset)
{
  if (node < 0) return;
  parse_tree[node].yr += offset;
  indent_yr (parse_tree, parse_tree[node].lchild, offset);
  indent_yr (parse_tree, parse_tree[node].rchild, offset);
}

void Pair_CFG_state_typing::indent_yl (Pair_CFG_parse_tree& parse_tree, int node, int offset)
{
  if (node < 0) return;
  parse_tree[node].yl += offset;
  indent_yl (parse_tree, parse_tree[node].lchild, offset);
  indent_yl (parse_tree, parse_tree[node].rchild, offset);
}

Pair_CFG_parse_tree Pair_CFG_state_typing::parse (const vector<int>& state_path) const
{
  return parse (state_path, 0, 0);
}

Pair_CFG_parse_tree Pair_CFG_state_typing::parse (const Pair_CFG_local_path& local_path) const
{
  // get parse tree
  const Pair_CFG_parse_tree parse_tree = parse (local_path.path, local_path.xstart, local_path.ystart);
  // check co-ords
  const Pair_CFG_branch& b = parse_tree.front();
  if (b.xl != local_path.xstart || b.xr != local_path.xstart + local_path.xlen
      || b.yl != local_path.ystart || b.yr != local_path.ystart + local_path.ylen) {
    CLOGERR << "State types:\n";
    for (int i = 0; i < states(); ++i)
      CL << " " << i << ": " << state_type_string (state_type[i]) << "\n";
    CLOGERR << "Parse tree:\n";
    parse_tree.show (CL);
    CLOGERR << "Local path:\n";
    local_path.show (CL);
    THROWEXPR ("Parse tree coords (xl=" << b.xl << ", xr=" << b.xr << ", yl=" << b.yl << ", yr=" << b.yr << ") don't match local path coords (xstart=" << local_path.xstart << ", xlen=" << local_path.xlen << ", ystart=" << local_path.ystart << ", ylen=" << local_path.ylen << ")");
  }
  // return
  return parse_tree;
}

Pair_CFG_parse_tree Pair_CFG_state_typing::parse (const vector<int>& state_path, int xstart, int ystart) const
{
  // during parsing, xr & yr vars in Pair_CFG_branch are measured as negative offsets from length
  // xl,xr,yl,yr keep track of the no. of residues emitted to left & right in X & Y for current branch only
  int xl = 0;
  int xr = 0;
  int yl = 0;
  int yr = 0;
  vector<int> xlen (1, (int) 0);  // for each branch: X-residues emitted by branch & all its children
  vector<int> ylen (1, (int) 0);  // for each branch: Y-residues emitted by branch & all its children
  Pair_CFG_parse_tree parse_tree (0, 0);
  // main parser loop
  for (int i = 0; i < (int) state_path.size(); ++i)
    {
      const int s = state_path[i];
      Pair_CFG_branch& current = parse_tree.back();
      current.path.push_back (s);
      if (i == 0)
	{
	  if (s != HMM_state_enum::Start) THROWEXPR ("No Start state");
	  continue;
	}
      if (s == HMM_state_enum::Start) THROWEXPR ("Unexpected Start state");
      if (s == HMM_state_enum::End)
	{
	  // find last parent with no right child
	  int c = parse_tree.size() - 1;  // child index
	  int p = -1;
	  while (c != 0)
	    {
	      p = parse_tree[c].parent;
	      xlen[p] += xlen[c];
	      ylen[p] += ylen[c];
	      if (parse_tree[p].rchild < 0) break;
	      // we just completed a right child, so update xr/yr for left-child (or xr/yl if type=BifurcRevY)
	      const int lc = parse_tree[p].lchild;
	      indent_xr (parse_tree, lc, xlen[c]);
	      if (parse_tree[p].revy)
		indent_yl (parse_tree, lc, ylen[c]);  // BifurcRevY
	      else
		indent_yr (parse_tree, lc, ylen[c]);  // Bifurc
	      c = p;
	    }
	  if (c > 0)
	    {
	      if (i == (int) state_path.size() - 1) THROWEXPR ("Nonempty stack");
	      // create right child
	      const Pair_CFG_branch& left = parse_tree[c];
	      if (parse_tree[p].revy)
		parse_tree.new_rchild (left.xl + xlen[c], left.xr, left.yl, left.yr + ylen[c], p);  // BifurcRevY
	      else
		parse_tree.new_rchild (left.xl + xlen[c], left.xr, left.yl + ylen[c], left.yr, p);  // Bifurc
	      xlen.push_back (0);
	      ylen.push_back (0);
	    }
	  else
	    if (i < (int) state_path.size() - 1) THROWEXPR ("Empty stack");

	  xl = xr = yl = yr = 0;
	}
      else  // s != End
	{
	  if (i == (int) state_path.size() - 1) THROWEXPR ("Unterminated branch");
	  const State_type t = state_type[s];
	  if (is_bifurc_type (t))
	    {
	      const int p = parse_tree.size() - 1;
	      parse_tree.new_lchild (current.xl + xl, current.xr + xr, current.yl + yl, current.yr + yr, p);
	      parse_tree[p].revy = t == BifurcRevY;
	      xlen.push_back (0);
	      ylen.push_back (0);
	      xl = xr = yl = yr = 0;
	    }
	  else  // emit state
	    {
	      if (t & EmitXL) { ++xl; ++xlen.back(); }
	      if (t & EmitXR) { ++xr; ++xlen.back(); }
	      if (t & EmitYL) { ++yl; ++ylen.back(); }
	      if (t & EmitYR) { ++yr; ++ylen.back(); }
	    }
	}
    }
  // calculate true xr & yr coords by adding subtracting from length
  // also add (xstart,ystart) to (xl,xr,yl,yr) so local alignments have correct co-ords
  const int total_xlen = xlen[0];
  const int total_ylen = ylen[0];

  if (CTAGGING(-1,PAIR_CFG_PARSE_TREE))
    CL << "In Pair_CFG_state_typing::parse: xstart=" << xstart << " ystart=" << ystart << " xlen=(" << xlen << ") ylen=(" << ylen << ")\n";
  for_contents (Pair_CFG_parse_tree, parse_tree, b)
    {
      if (CTAGGING(-1,PAIR_CFG_PARSE_TREE))
	CL << "Offsetting Pair_CFG_branch coords: xl=" << b->xl << " xr=" << b->xr << " yl=" << b->yl << " yr=" << b->yr << "\n";
      b->xl += xstart;
      b->yl += ystart;
      b->xr = xstart + total_xlen - b->xr;
      b->yr = ystart + total_ylen - b->yr;
    }
  return parse_tree;
}

void Pair_CFG_state_typing::reset_state_types()
{
  state_type = vector<State_type> (states(), (State_type) Undefined);
  bifurc = vector<Bifurcation> (states(), Bifurcation());
}

bool Pair_CFG_state_typing::same_state_types (const Pair_CFG_state_typing& t) const
{
  return state_type == t.state_type && bifurc == t.bifurc;
}

vector<vector<Pair_CFG_state_type_enum::Bifurcation_left_parent> > Pair_CFG_state_typing::left_parent() const
{
  vector<vector<Bifurcation_left_parent> > lp (states());
  for (int s = 0; s < states(); ++s)
    if (is_bifurc_type (state_type[s]))
      lp[bifurc[s].r].push_back (Bifurcation_left_parent (bifurc[s].l, s));
  return lp;
}

vector<vector<Pair_CFG_state_type_enum::Bifurcation_right_parent> > Pair_CFG_state_typing::right_parent() const
{
  vector<vector<Bifurcation_right_parent> > rp (states());
  for (int s = 0; s < states(); ++s)
    if (is_bifurc_type (state_type[s]))
      rp[bifurc[s].l].push_back (Bifurcation_right_parent (bifurc[s].r, s));
  return rp;
}

vector<int> Pair_CFG_state_typing::unparse (const Pair_CFG_parse_tree& parse_tree) const
{
  vector<int> state_path;
  for_const_contents (Pair_CFG_parse_tree, parse_tree, b)
    state_path.insert (state_path.end(), b->path.begin(), b->path.end());
  return state_path;
}

bool Pair_CFG_state_typing::test_state_types_valid() const
{
  if (state_type.size() != bifurc.size()) return 0;
  for (int s = 0; s < states(); ++s)
    {
      const State_type& t = state_type[s];
      if (!test_state_type_defined(t))
	return 0;
      const Bifurcation& b = bifurc[s];
      if (is_bifurc_type(t))
	{
	  if (b.l < 0 || b.l >= states() || b.r < 0 || b.r >= states()) return 0;
	  if (is_emit_type (state_type[b.l]) || is_emit_type (state_type[b.r])) return 0;
	}
      else
	if (!b.null()) return 0;
    }
  return 1;
}

bool Pair_CFG_state_typing::test_branch_coords_consistent (const Pair_CFG_parse_tree& parse_tree) const
{
  if (!parse_tree.test_connections()) return 0;
  for_const_contents (Pair_CFG_parse_tree, parse_tree, b)
    {
      // find number of residues emitted in each direction
      int xl = 0;
      int xr = 0;
      int yl = 0;
      int yr = 0;
      for (int i = 0; i < (int) b->path.size(); ++i)
	if (b->path[i] >= 0)
	  {
	    const State_type t = state_type[b->path[i]];
	    if (is_bifurc_type (t) && i != (int) b->path.size() - 1) return 0;
	    if (t & EmitXL) ++xl;
	    if (t & EmitXR) ++xr;
	    if (t & EmitYL) ++yl;
	    if (t & EmitYR) ++yr;
	  }
      // test consistency
      const int s = b->path.back();
      if (s == HMM_state_enum::End)
	{
	  if (b->has_children()) return 0;
	  if (b->xl + xl != b->xr - xr || b->yl + yl != b->yr - yr) return 0;
	}
      else
	{
	  const State_type t = state_type[s];
	  if (!is_bifurc_type (t)) return 0;  // each branch must end with an End or a bifurcation
	  if (!b->has_children()) return 0;
	  const Pair_CFG_branch& lchild = parse_tree[b->lchild];
	  const Pair_CFG_branch& rchild = parse_tree[b->rchild];
	  if (lchild.path[0] != bifurc[s].l || rchild.path[0] != bifurc[s].r) return 0;
	  if (lchild.xl != b->xl + xl || rchild.xr != b->xr - xr) return 0;
	  if (t != (b->revy ? BifurcRevY : Bifurc)) return 0;
	  if (b->revy
	      ? (rchild.yl != b->yl + yl || lchild.yr != b->yr - yr)
	      : (lchild.yl != b->yl + yl || rchild.yr != b->yr - yr))
	    return 0;
	}
    }
  return 1;
}

bool Pair_CFG_state_typing::is_single_CFG() const
{
  for (int state = 0; state < states(); ++state)
    {
      const int type = state_type[state];
      if (emit_dyl(type) || emit_dyr(type)) return FALSE;
    }
  return TRUE;
}

void Pair_CFG_scores::init_emit (int state, State_type type, Score sc)
{
  Pair_CFG<Score>::init_emit (state, type, sc);
}

vector<int> Pair_CFG_scores::emit_states() const
{
  vector<int> result;
  for (int s = 0; s < states(); ++s) if (is_emit_type (state_type[s])) result.push_back(s);
  return result;
}

vector<int> Pair_CFG_scores::nonemit_states_unsorted() const
{
  vector<int> result;
  for (int s = 0; s < states(); ++s) if (!is_emit_type (state_type[s])) result.push_back(s);
  return result;
}

vector<int> Pair_CFG_scores::nonemit_states() const
{
  // the following #if'd-out line ensures bifurcation states are properly sorted (i.e. filled after their target states)
  //  if bifurcations involving empty subsequences are allowed.
  // however, this requires that bifurcation states never bifurcate back to themselves
  //  (since bifurcations are considered to be null transitions).
  // practically, it's more useful to have self-looping bifurcations than empty-subsequence bifurcations,
  //  so this code is not compiled in, but it's left here in case, one day, I decide to make this optional.
  //
  // Addendum (29/1/2004): ran into problems with the TKF structure tree grammar, which *does* have self-looping empty bifurcations.
  // Fixed this by patching the grammar to add direct transitions that reproduce the effect of these bifurcations.
#if FALSE
  ((Pair_CFG_scores*) this) -> add_fake_bifurcation_transitions();  // a hacky cast for a hacky technique
#endif /* FALSE */
  const vector<int> unsorted = nonemit_states_unsorted();
  return Transition_methods::topological_sort (*this, unsorted);
}

void Pair_CFG_scores::zero_outgoing_bifurcation_transitions()
{
  for (int s = 0; s < states(); ++s)
    if (is_bifurc_type (state_type[s]))
      for (int d = 0; d < states(); ++d)
	transition (s, d) = -InfinityScore;
}

void Pair_CFG_scores::add_fake_bifurcation_transitions()
{
  zero_outgoing_bifurcation_transitions();
  for (int s = 0; s < states(); ++s)
    if (is_bifurc_type (state_type[s]))
      {
	transition (s, bifurc[s].l) = 0;
	transition (s, bifurc[s].r) = 0;
      }
}

vector<vector<int> > Pair_CFG_scores::incoming_states() const
{
  ((Pair_CFG_scores*) this) -> zero_outgoing_bifurcation_transitions();  // cast away const.. hacky
  return Transition_methods::incoming_states (*this);
}

vector<vector<int> > Pair_CFG_scores::selected_outgoing_states (const vector<int>& selection) const
{
  ((Pair_CFG_scores*) this) -> zero_outgoing_bifurcation_transitions();  // cast away const.. hacky
  return Transition_methods::selected_outgoing_states (*this, selection);
}

vector<vector<int> > Pair_CFG_scores::selected_incoming_states (const vector<int>& selection) const
{
  ((Pair_CFG_scores*) this) -> zero_outgoing_bifurcation_transitions();  // cast away const.. hacky
  return Transition_methods::selected_incoming_states (*this, selection);
}

Score Pair_CFG_scores::path_transition_score (const Pair_CFG_parse_tree& parse_tree) const
{
  Score sc = 0;
  for_const_contents (Pair_CFG_parse_tree, parse_tree, b)
    ScorePMulAcc (sc, Transition_methods::path_transition_score (*this, b->path));
  return sc;
}

Score Pair_CFG_scores::path_emit_score (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy) const
{
  Score sc = 0;
  for_const_contents (Pair_CFG_parse_tree, parse_tree, b)
    {
      int xl = b->xl;
      int xr = b->xr;
      int yl = b->yl;
      int yr = b->yr;
      for_const_contents (vector<int>, b->path, s)
	if (*s >= 0)
	  {
	    const State_type t = state_type[*s];
	    if (is_emit_type(t))
	      {
		int emit_idx = 0;
		if (t & EmitXL) emit_idx += emit_xl_mul(t) * dsqx[xl++];
		if (t & EmitXR) emit_idx += emit_xr_mul(t) * dsqx[--xr];
		if (t & EmitYL) emit_idx += emit_yl_mul(t) * dsqy[yl++];
		if (t & EmitYR) emit_idx += emit_yr_mul(t) * dsqy[--yr];
		const Score emit_sc = emit[*s][emit_idx];
		ScorePMulAcc (sc, emit_sc);
	      }
	  }
    }
  return sc;
}

Score Pair_CFG_scores::path_score (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy) const
{
  return ScorePMul (path_transition_score (parse_tree), path_emit_score (parse_tree, dsqx, dsqy));
}

Pair_CFG_scores::XYMask Pair_CFG_scores::get_mask (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy, const set<int>& mask_states) const
{
  XYMask mask;
  mask.xmask = vector<bool> (dsqx.size());
  mask.ymask = vector<bool> (dsqy.size());
  for_const_contents (Pair_CFG_parse_tree, parse_tree, b)
    {
      int xl = b->xl;
      int xr = b->xr;
      int yl = b->yl;
      int yr = b->yr;
      for_const_contents (vector<int>, b->path, s)
	if (*s >= 0)
	  {
	    const bool masked = mask_states.find(*s) != mask_states.end();
	    const State_type t = state_type[*s];
	    if (is_emit_type(t))
	      {
		if (t & EmitXL) mask.xmask[xl++] = masked;
		if (t & EmitXR) mask.xmask[--xr] = masked;
		if (t & EmitYL) mask.ymask[yl++] = masked;
		if (t & EmitYR) mask.ymask[--yr] = masked;
	      }
	  }
    }
  return mask;
}

GFF Pair_CFG_scores::ymask2gff (const XYMask& mask, const Named_profile& xprof, const Named_profile& yprof) const
{
  int xstart, xend, ystart, yend;
  for (xstart = 0; xstart < (int) mask.xmask.size(); ++xstart) if (mask.xmask[xstart]) break;
  for (xend = (int) mask.xmask.size() - 1; xend >= 0; --xend) if (mask.xmask[xend]) break;
  for (ystart = 0; ystart < (int) mask.ymask.size(); ++ystart) if (mask.ymask[ystart]) break;
  for (yend = (int) mask.ymask.size() - 1; yend >= 0; --yend) if (mask.ymask[yend]) break;
  GFF gff;
  gff.seqname = yprof.name;
  gff.start = ystart + 1;
  gff.end = yend + 1;
  gff.feature << xprof.name << "/" << xstart+1 << "-" << xend+1 << ":";
  for (int x = xstart; x <= xend; ++x) gff.feature << xprof.seq[x];
  for (int y = ystart; y <= yend; ++y) gff.group << yprof.seq[y];
  return gff;
}

void Pair_CFG_scores::show_score_breakdown (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy, ostream& out) const
{
  typedef pair<int,int> Branch_indent;
  stack<Branch_indent> next_branch;
  next_branch.push (Branch_indent (0, 0));
  while (!next_branch.empty())
    {
      const Branch_indent bi = next_branch.top();
      const Pair_CFG_branch& b = parse_tree[bi.first];
      const int indent = bi.second;
      next_branch.pop();
      int xl = b.xl;
      int xr = b.xr;
      int yl = b.yl;
      int yr = b.yr;
      for (int i = 0; i < indent; ++i)
	out << ' ';
      out << "[" << b.xl << ".." << b.xr << "," << b.yl << ".." << b.yr << "]";
      for (int pos = 0; pos < (int) b.path.size(); ++pos)
	{
	  const int d = b.path[pos];
	  if (pos > 0)
	    {
	      const int s = b.path[pos-1];
	      const Score trans_sc = transition (s, d);
	      out << ' ';
	      if (trans_sc >= 0) out << '+';
	      ShowScore(trans_sc,out);
	      out << '(' << s << "=>" << d << ')';
	    }
	  if (d >= 0)
	    {
	      const State_type t = state_type[d];
	      if (is_emit_type(t))
		{
		  int emit_idx = 0;
		  vector<int> emitted (4, -1);
		  if (t & EmitXL) emit_idx += emit_xl_mul(t) * (emitted[0] = dsqx[xl++]);
		  if (t & EmitXR) emit_idx += emit_xr_mul(t) * (emitted[1] = dsqx[--xr]);
		  if (t & EmitYL) emit_idx += emit_yl_mul(t) * (emitted[2] = dsqy[yl++]);
		  if (t & EmitYR) emit_idx += emit_yr_mul(t) * (emitted[3] = dsqy[--yr]);
		  sstring emit_string ("....");
		  for (int i = 0; i < 4; ++i)
		    if (emitted[i] >= 0)
		      emit_string[i] = alphabet().int2char_uc (emitted[i]);
		  const Score emit_sc = emit[d][emit_idx];
		  out << ' ';
		  if (emit_sc >= 0) out << '+';
		  ShowScore(emit_sc,out);
		  out << '(' << emit_string << ')';
		}
	    }
	}
      out << "\n";
      if (b.has_children())
	{
	  next_branch.push (Branch_indent (b.rchild, indent+1));
	  next_branch.push (Branch_indent (b.lchild, indent+1));
	}
    }
}

void Pair_CFG_counts::add_counts_from_parse_tree (const Pair_CFG_scores& cfg, const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy)
{
  log_likelihood += Score2Nats (cfg.path_score (parse_tree, dsqx, dsqy));

  for_const_contents (Pair_CFG_parse_tree, parse_tree, b)
    {
      Transition_methods::add_transition_counts_from_path (*this, b->path);

      int xl = b->xl;
      int xr = b->xr;
      int yl = b->yl;
      int yr = b->yr;
      for_const_contents (vector<int>, b->path, s)
	if (*s >= 0)
	  {
	    const int t = state_type[*s];
	    if (is_emit_type(t))
	      {
		int emit_idx = 0;
		if (t & EmitXL) emit_idx += emit_xl_mul(t) * dsqx[xl++];
		if (t & EmitXR) emit_idx += emit_xr_mul(t) * dsqx[--xr];
		if (t & EmitYL) emit_idx += emit_yl_mul(t) * dsqy[yl++];
		if (t & EmitYR) emit_idx += emit_yr_mul(t) * dsqy[--yr];
		emit[*s][emit_idx] += 1.0;
	      }
	  }
    }
}

set<int> Pair_CFG_branch::get_paired_states (const vector<State_type>& state_type)
{
  set<int> paired_states;
  for (int s = 0; s < (int) state_type.size(); ++s)
    {
      const int t = state_type[s];
      if ((t & EmitXL && t & EmitXR) || (t & EmitYL && t & EmitYR))
	paired_states.insert (s);
    }
  return paired_states;
}

Pair_CFG_alignment Pair_CFG_branch::alignment (const vector<Pair_CFG_branch>& parse_tree, const vector<State_type>& state_type, const Named_profile& npx, const Named_profile& npy) const
{
  const set<int> paired_states = get_paired_states (state_type);
  return alignment (parse_tree, state_type, paired_states, npx, npy);
}

Pair_CFG_alignment Pair_CFG_branch::alignment (const vector<Pair_CFG_branch>& parse_tree, const vector<State_type>& state_type, const set<int>& paired_states, const Named_profile& npx, const Named_profile& npy) const
{
  sstring xlfold, ylfold, xrfold, yrfold;
  vector<State_type> parse;
  int xlpos = xl;
  int xrpos = xr;
  int ylpos = yl;
  int yrpos = yr;
  for (int i = 0; i < (int) path.size(); ++i)
    {
      const int s = path[i];
      if (s < 0) continue;
      const State_type t = state_type[s];
      const int emitl = t & (EmitXL | EmitYL);
      const int emitr = t & (EmitXR | EmitYR);
      int xlc = 0, xrc = 0, ylc = 0, yrc = 0;
      if (t & EmitXL) xlc = npx.dsq[xlpos++];
      if (t & EmitXR) xrc = npx.dsq[--xrpos];
      if (t & EmitYL) ylc = npy.dsq[ylpos++];
      if (t & EmitYR) yrc = npy.dsq[--yrpos];
      const bool is_paired_state = paired_states.find(s) != paired_states.end();
      const bool is_xpaired = is_paired_state && (t & EmitXL) && (t & EmitXR);
      const bool is_ypaired = is_paired_state && (t & EmitYL) && (t & EmitYR);
      const int wcx = is_xpaired && xlc == CFG_alphabet.complement (xrc);
      const int wcy = is_ypaired && ylc == CFG_alphabet.complement (yrc);
      if (emitl)
	{
	  xlfold << (t & EmitXL ? (is_xpaired ? (wcx ? (char) FOLD_LCHAR : (char) NONWC_LCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  ylfold << (t & EmitYL ? (is_ypaired ? (wcy ? (char) FOLD_LCHAR : (char) NONWC_LCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	}
      if (emitr)
	{
	  xrfold << (t & EmitXR ? (is_xpaired ? (wcx ? (char) FOLD_RCHAR : (char) NONWC_RCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  yrfold << (t & EmitYR ? (is_ypaired ? (wcy ? (char) FOLD_RCHAR : (char) NONWC_RCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	}
      if (t != Null)
	parse.push_back (t);
    }
  if (has_children())
    {
      const Pair_CFG_alignment lchild_align = parse_tree[lchild].alignment (parse_tree, state_type, paired_states, npx, npy);
      const Pair_CFG_alignment rchild_align = parse_tree[rchild].alignment (parse_tree, state_type, paired_states, npx, npy);
      xlfold << lchild_align.xfold << rchild_align.xfold;
      ylfold << lchild_align.yfold << rchild_align.yfold;
      parse.insert (parse.end(), lchild_align.parse.begin(), lchild_align.parse.end());
      parse.insert (parse.end(), rchild_align.parse.begin(), rchild_align.parse.end());
      if (revy) THROWEXPR ("BifurcRevY state unimplemented for alignment views");
    }
  else
    parse.push_back (Null);  // Null in parse indicates a leaf node
  reverse (xrfold.begin(), xrfold.end());
  reverse (yrfold.begin(), yrfold.end());
  xlfold << xrfold;
  ylfold << yrfold;
  Pair_CFG_alignment align (npx, npy);
  align.parse.swap (parse);
  align.xfold.swap (xlfold);
  align.yfold.swap (ylfold);
  align.xstart = xl;
  align.ystart = yl;
  return align;
}

Pair_CFG_alignment Pair_CFG_parse_tree::alignment (const vector<State_type>& state_type, const set<int>& paired_states, const Named_profile& npx, const Named_profile& npy) const
{
  if (size() == 0)
    {
      Pair_CFG_alignment align (npx, npy);
      return align;
    }
  return front().alignment (*this, state_type, paired_states, npx, npy);
}

Pair_CFG_alignment Pair_CFG_parse_tree::alignment (const vector<State_type>& state_type, const Named_profile& npx, const Named_profile& npy) const
{
  const set<int> paired_states = Pair_CFG_branch::get_paired_states (state_type);
  return alignment (state_type, paired_states, npx, npy);
}

Pair_CFG_alignment::Pair_CFG_alignment (const Named_profile& npx, const Named_profile& npy)
  : npx (npx), npy (npy)
{
  xseqlen = npx.size();
  yseqlen = npy.size();
}

sstring RNA_pair_path::xfold_ungapped() const
{
  sstring fold;
  for_const_contents (sstring, xfold, c)
    if (*c != GAP_CHAR)
      fold.push_back (*c);
  return fold;
}

sstring RNA_pair_path::yfold_ungapped() const
{
  sstring fold;
  for_const_contents (sstring, yfold, c)
    if (*c != GAP_CHAR)
      fold.push_back (*c);
  return fold;
}

int RNA_pair_path::xfold_ungapped_len() const
{
  int len = 0;
  for_const_contents (sstring, xfold, c)
    if (*c != GAP_CHAR)
      ++len;
  return len;
}

int RNA_pair_path::yfold_ungapped_len() const
{
  int len = 0;
  for_const_contents (sstring, yfold, c)
    if (*c != GAP_CHAR)
      ++len;
  return len;
}

Pairwise_path RNA_pair_path::pairwise_path (bool global_pad) const
{
  if (xfold.size() != yfold.size())
    THROWEXPR ("Gapped fold strings have different lengths: X='" << xfold << "', Y='" << yfold << "'");

  Pairwise_path path;
  if (global_pad)
    {
      for (int i = 0; i < xstart; ++i) path.append_column (1, 0);
      for (int i = 0; i < ystart; ++i) path.append_column (0, 1);
    }
  for (int i = 0; i < (int) xfold.size(); ++i)
    path.append_column (xfold[i] != GAP_CHAR, yfold[i] != GAP_CHAR);
  if (global_pad)
    {
      for (int i = xstart + xfold_ungapped_len(); i < xseqlen; ++i) path.append_column (1, 0);
      for (int i = ystart + yfold_ungapped_len(); i < yseqlen; ++i) path.append_column (0, 1);
    }
  return path;
}

void RNA_pair_path::add_pairwise_path (Pair_envelope& env, bool allow) const
{
  const Pairwise_path path = pairwise_path (TRUE);
  env.add_pairwise_path (path, allow);
}

Local_fold_string RNA_pair_path::local_xfold() const
{
  return Local_fold_string (xfold_ungapped(), xstart, xseqlen);
}

Local_fold_string RNA_pair_path::local_yfold() const
{
  return Local_fold_string (yfold_ungapped(), ystart, yseqlen);
}

Stockade Pair_CFG_alignment::stockade (bool use_NSE) const
{
  // get pairwise path
  const Pairwise_path path = pairwise_path (FALSE);
  const int xlen = path.count_steps_in_row(0);
  const int ylen = path.count_steps_in_row(1);
  // create Stockade
  Stockade stock (2);
  // initialise Named_profiles
  stock.np[0] = npx.subseq (xstart, xlen);
  stock.np[1] = npy.subseq (ystart, ylen);
  // update names in Stockade
  stock.align.row_name[0] = stock.np[0].name;
  stock.align.row_name[1] = stock.np[1].name;
  if (use_NSE)
    {
      stock.align.row_name[0] << '/' << xstart << '-' << xstart+xlen-1;	
      stock.align.row_name[1] << '/' << ystart << '-' << ystart+ylen-1;
    }
  stock.align.update_row_index();
  // set alignment path
  stock.align.path = path;
  // set fold strings
  stock.align.set_gr_annot (stock.align.row_name[0], Stockholm_secondary_structure_tag, xfold);
  stock.align.set_gr_annot (stock.align.row_name[1], Stockholm_secondary_structure_tag, yfold);
  // return Stockade
  return stock;
}

void Pair_CFG_alignment::show (ostream& o) const
{
  save_flags (o);
  left_align (o);

  if (xfold.size() != yfold.size()) THROWEXPR ("Alignment X- and Y-fold lengths don't match");

  sstring gappedx (xfold.size());
  sstring gappedy (yfold.size());
  int i, xpos, ypos;

  for (i = 0, xpos = xstart; i < (int) xfold.size(); ++i)
    gappedx[i] = xfold[i] == GAP_CHAR ? GAP_CHAR : npx.seq[xpos++];

  for (i = 0, ypos = ystart; i < (int) yfold.size(); ++i)
    gappedy[i] = yfold[i] == GAP_CHAR ? GAP_CHAR : npy.seq[ypos++];

  sstring ss_cons (xfold.size());
  sstring ps_cons (xfold.size());
  for (i = 0; i < (int) xfold.size(); ++i)
    {
      ps_cons[i] = gappedx[i] == gappedy[i] ? gappedx[i] : '.';
      ss_cons[i] = is_lchar(xfold[i]) && is_lchar(yfold[i]) ? FOLD_LCHAR : (is_rchar(xfold[i]) && is_rchar(yfold[i]) ? FOLD_RCHAR : NOFOLD_CHAR);
    }

  sstring xname, yname, bits, xss, yss, ss_cons_name, ps_cons_name;
  xname << npx.name << '/' << xstart+1 << '-' << xpos;
  yname << npy.name << '/' << ystart+1 << '-' << ypos;
  bits << Stockholm_file_annotation << ' ' << Stockholm_bit_score_tag;
  xss << Stockholm_sequence_column_annotation << ' ' << xname << ' ' << Stockholm_secondary_structure_tag;
  yss << Stockholm_sequence_column_annotation << ' ' << yname << ' ' << Stockholm_secondary_structure_tag;
  ps_cons_name << Stockholm_column_annotation << ' ' << StockholmConsensus(Stockholm_primary_sequence_tag);
  ss_cons_name << Stockholm_column_annotation << ' ' << StockholmConsensus(Stockholm_secondary_structure_tag);
  const int name_len = max (max (xss.size(), yss.size()), max (ss_cons_name.size(), bits.size()));

  o << "# " << Stockholm_format_identifier << " " << Stockholm_version_identifier << "\n";

  o.width (name_len);
  o << xss << " " << xfold << "\n";
  o.width (name_len);
  o << xname << " " << gappedx << "\n";

  o.width (name_len);
  o << ps_cons_name << " " << ps_cons << "\n";

  o.width (name_len);
  o << yname << " " << gappedy << "\n";
  o.width (name_len);
  o << yss << " " << yfold << "\n";

  o.width (name_len);
  o << ss_cons_name << " " << ss_cons << "\n";

  o.width (name_len);
  o << bits << " " << Score2Bits (score) << "\n";

  o << Stockholm_alignment_separator << "\n";

  restore_flags (o);
}


// Elimination of null bifurcations (transformation A->B) and null state paths (transformation B->C).

// returns null Inside probabilities, P(empty_sequence|state), indexed by state (excluding Start & End)
// currently computes these probabilities via a simple iterative lower-bound, heuristically set at 100*N*log(N) iterations, where N is the number of null states
#define EMPTY_PROBS_HEURISTIC 100
vector<Prob> Pair_CFG_scores::empty_probs() const
{
  vector<Prob> p_empty (states(), 0.);
  const vector<int> null = nonemit_states_unsorted();
  const int max_iter = (int) (((double) null.size()) * Math_fn::math_log ((double) null.size()) * EMPTY_PROBS_HEURISTIC);   // completely arbitrary heuristic formula for number of iterations...
  CTAG(4,CFG_ELIM) << "Using iterative heuristic to find null subtree probabilities for Pair SCFG (" << max_iter << " iterations)\n";
  for (int iter = 0; iter < max_iter; ++iter)
    for_const_contents (vector<int>, null, s)
      {
	const State_type t = state_type[*s];
	if (is_bifurc_type (t))
	  p_empty[*s] = p_empty[bifurc[*s].l] * p_empty[bifurc[*s].r];
	else
	  {
	    Prob p = Score2Prob (transition (*s, End));
	    for_const_contents (vector<int>, null, d)
	      p += Score2Prob (transition (*s, *d)) * p_empty[*d];
	    p_empty[*s] = p;
	  }
      }
  CTAG(4,CFG_ELIM) << "Finished " << max_iter << " iterations; p_empty = (" << p_empty << ")\n";
  return p_empty;
}

vector<int> Pair_CFG_scores::sample_null_subtree (const vector<Prob>& p_empty,  // a co-ords
						  const Pair_CFG_scores& a_cfg,
						  int root_state)  // a co-ords
{
  // currently implemented recursively
  vector<int> subtree (1, (int) root_state);
  const Pair_CFG_state_type_enum::State_type t = a_cfg.state_type[root_state];
  if (a_cfg.is_bifurc_type (t))
    {
      const vector<int>& l_subtree = sample_null_subtree (p_empty, a_cfg, a_cfg.bifurc[root_state].l);
      const vector<int>& r_subtree = sample_null_subtree (p_empty, a_cfg, a_cfg.bifurc[root_state].r);
      subtree.insert (subtree.end(), l_subtree.begin(), l_subtree.end());
      subtree.insert (subtree.end(), r_subtree.begin(), r_subtree.end());
    }
  else
    {
      vector<int> dest = a_cfg.nonemit_states_unsorted();
      dest.push_back (Grammar_state_enum::End);

      while (root_state != Grammar_state_enum::End)
	{
	  vector<Prob> dest_p;
	  for_const_contents (vector<int>, dest, d)
	    dest_p.push_back (Score2Prob (a_cfg.transition (root_state, *d)) * (*d == Grammar_state_enum::End ? 1. : p_empty[*d]));
	  root_state = dest[Rnd::choose (dest_p)];
	  subtree.push_back (root_state);
	}
    }
  return subtree;
}

Pair_CFG_scores Pair_CFG_scores::eliminate_null_states_and_bifurcations (vector<Prob>& empty_probs,  // a co-ordinates
									 Concrete_transition_probs& tp_orig,  // b co-ordinates
									 Concrete_transition_probs& tp_elim,  // b co-ordinates
									 vector<int>& b2a_state,  // b->a mapping
									 vector<int>& c2b_state) const
{
  vector<int> bifurc_prec_states;
  const Pair_CFG_scores b_cfg = add_null_bifurc_transitions (empty_probs, b2a_state, bifurc_prec_states);
  const Pair_CFG_scores c_cfg = b_cfg.eliminate_null_states (tp_orig, tp_elim, bifurc_prec_states, c2b_state);
  return c_cfg;
}

vector<int> Pair_CFG_scores::sample_eliminated_states_and_bifurcations (const vector<Prob>& empty_probs,  // a co-ords
									const Transition_probs& tp_orig,  // b co-ords
									const Transition_probs& tp_elim,  // b co-ords
									const vector<int>& b2a_state,  // b->a
									const vector<int>& c2b_state,  // c->b
									const Pair_CFG_scores& a_cfg,
									const Pair_CFG_scores& b_cfg,
									const vector<int>& c_path)
{
  const vector<int> b_path = sample_eliminated_states (tp_orig, tp_elim, c2b_state, b_cfg, c_path);
  const vector<int> a_path = sample_eliminated_bifurcations (empty_probs, b2a_state, a_cfg, b_cfg, b_path);
  return a_path;
}

#define CFG_NULBIF_SUFFIX ".nulbif"
Pair_CFG_scores Pair_CFG_scores::add_null_bifurc_transitions (vector<Prob>& p_empty,  // a co-ordinates
							      vector<int>& b2a_state,  // b->a
							      vector<int>& bifurc_prec_states) const  // b co-ordinates
{
  // create a new null state preceding every bifurcation state; move bifurcation states to the end
  // add transitions from null state representing "empty" bifurcations

  // get null subtree probs
  p_empty = empty_probs();

  // create a<-->b state maps (a-->b contains bifurc states only)
  b2a_state.clear();
  bifurc_prec_states.clear();
  map<int,int> a2b_state;
  for (int s = 0; s < states(); ++s)
    b2a_state.push_back (s);
  for (int s = 0; s < states(); ++s)
    if (is_bifurc_type (state_type[s]))
      {
	a2b_state[s] = b2a_state.size();
	b2a_state.push_back (s);
	bifurc_prec_states.push_back (s);
      }

  // create pair SCFG b
  Pair_CFG_scores b_cfg (b2a_state.size());
  b_cfg.name.clear();
  b_cfg.name << name << CFG_NULBIF_SUFFIX;

  // set up b
  for (int s = 0; s < states(); ++s)
    {
      b_cfg.state_name[s] = state_name[s];
      if (is_bifurc_type (state_type[s]))
	{
	  // bifurcation state, so get new state index
	  const int s2 = a2b_state[s];
	  // set state name
	  b_cfg.state_name[s2].clear();
	  b_cfg.state_name[s2] << state_name[s] << CFG_NULBIF_SUFFIX;
	  // set state types
	  b_cfg.state_type[s] = Null;
	  b_cfg.state_type[s2] = state_type[s];
	  // copy bifurcation info
	  b_cfg.bifurc[s2] = bifurc[s];
	  // add transitions from null state representing empty bifurcations
	  const Score l_empty_sc = Prob2Score (p_empty[bifurc[s].l]);
	  const Score r_empty_sc = Prob2Score (p_empty[bifurc[s].r]);
	  ScorePSumAcc (b_cfg.transition (s, bifurc[s].l), r_empty_sc);
	  ScorePSumAcc (b_cfg.transition (s, bifurc[s].r), l_empty_sc);
	  // add transition from null state to end state, representing double-empty bifurcation
	  ScorePSumAcc (b_cfg.transition (s, End), ScorePMul (l_empty_sc, r_empty_sc));
	  // add probability-one transition from null state to actual bifurcation state for non-empty bifurcations
	  // (caller must call Fold_envelope::remove_empty_bifurcations() to enforce non-empty constraint on DP algorithms when Pair SCFGs with eliminated null bifurcations are used)
	  b_cfg.transition (s, s2) = 0;
	  // copy start transition
	  b_cfg.start[s] = start[s];
	}
      else
	{
	  // copy everything over
	  b_cfg.state_name[s] = state_name[s];
	  b_cfg.state_type[s] = state_type[s];
	  b_cfg.emit[s] = emit[s];
	  b_cfg.start[s] = start[s];
	  b_cfg.end[s] = end[s];
	  for (int d = 0; d < states(); ++d)
	    b_cfg.transition(s,d) = transition(s,d);
	  // copy silly meta stuff too
	  b_cfg.xlmeta_idx[s] = xlmeta_idx[s];
	  b_cfg.xrmeta_idx[s] = xrmeta_idx[s];
	  b_cfg.ylmeta_idx[s] = ylmeta_idx[s];
	  b_cfg.ylmeta_idx[s] = ylmeta_idx[s];
	}
    }

  // return
  return b_cfg;
}

// B->C: eliminates all null states, factoring paths through null states into direct transitions between emit/bifurcation states
// (null states that are bifurcation targets are preserved, but are left inaccessible except through those bifurcations)
// method call on b_cfg returns c_cfg
#define CFG_ELIM_SUFFIX ".elim"
Pair_CFG_scores Pair_CFG_scores::eliminate_null_states (Concrete_transition_probs& tp_orig,  // b co-ordinates
							Concrete_transition_probs& tp_elim,  // b co-ordinates
							const vector<int>& bifurc_prec_states,  // b co-ordinates
							vector<int>& c2b_state) const  // c->b mapping (state indices change due to removal of null states)
{
  // get emit & nonemit states
  vector<int> null_states, emit_and_bifurc_states;
  for (int s = 0; s < states(); ++s)
    if (state_type[s] == Null)
      null_states.push_back (s);
    else
      emit_and_bifurc_states.push_back (s);

  // get tp_orig
  tp_orig = Transition_methods::score2prob (*this);

  // This all gets very convoluted because we have to handle sum-over-paths transitions to End separately from transitions to emit
  //  (transitions from bifurc_prec_states to non-End states can only be used for nonempty sequence, i.e. not on paths to End).
  // Also, Transition_methods::eliminate returns transitions from null states in one matrix (loop_exit) and from emit states in another (tp_elim).
  // Thus, there are four matrices we need to merge into one.

  // sum over paths through null states
  Concrete_transition_probs loop_exit_nonempty (states(), 0.);

  CTAG(1,CFG_ELIM) << "In Pair_CFG_scores::eliminate_null_states: allowing non-empty transitions\n";
  Concrete_transition_probs tp_elim_nonempty = Transition_methods::eliminate (tp_orig, null_states, false, &loop_exit_nonempty);

  // sum over paths through null states, EXCLUDING nonempty transitions (i.e. from bifurc_prec states to everything except the End state)
  Concrete_transition_probs tp_orig_empty (tp_orig);
  for_const_contents (vector<int>, bifurc_prec_states, bp)
    for (int d = 0; d < states(); ++d)
      tp_orig_empty.transition (*bp, d) = 0.;

  Concrete_transition_probs loop_exit_empty (states(), 0.);

  CTAG(1,CFG_ELIM) << "In Pair_CFG_scores::eliminate_null_states: disallowing non-empty transitions\n";
  Concrete_transition_probs tp_elim_empty = Transition_methods::eliminate (tp_orig_empty, null_states, false, &loop_exit_empty);

  // set up tp_elim, taking transitions from the following places:
  //  Start->emit from tp_elim_nonempty
  //   emit->emit from tp_elim_nonempty
  //   null->null from tp_elim_nonempty
  //   emit->null from tp_elim_nonempty
  //   null->emit from loop_exit_nonempty
  //   emit->End  from tp_elim_empty
  //   null->End  from loop_exit_empty

  // Start->emit, emit->emit, null->null, emit->null: tp_elim_nonempty
  tp_elim = tp_elim_nonempty;

  // null->emit: loop_exit_nonempty
  for_const_contents (vector<int>, null_states, s)
    for_const_contents (vector<int>, emit_and_bifurc_states, d)
    tp_elim.transition (*s, *d) = loop_exit_nonempty.transition (*s, *d);

  // emit->End: tp_elim_empty
  for_const_contents (vector<int>, emit_and_bifurc_states, s)
    tp_elim.transition (*s, End) = tp_elim_empty.transition (*s, End);

  // null->End: loop_exit_empty
  for_const_contents (vector<int>, null_states, s)
    tp_elim.transition (*s, End) = loop_exit_empty.transition (*s, End);

  // log
  if (CTAGGING(1,CFG_ELIM))
    {
      CL << "null_states: " << null_states << '\n';
      CL << "emit_and_bifurc_states: " << emit_and_bifurc_states << '\n';
      CL << "bifurc_prec_states: " << bifurc_prec_states << '\n';
      CL << "tp_orig "; tp_orig.show_transitions (CL);
      CL << "tp_elim_nonempty "; tp_elim_nonempty.show_transitions (CL);
      CL << "loop_exit_nonempty "; loop_exit_nonempty.show_transitions (CL);
      CL << "tp_orig_empty "; tp_orig_empty.show_transitions (CL);
      CL << "tp_elim_empty "; tp_elim_empty.show_transitions (CL);
      CL << "loop_exit_empty "; loop_exit_empty.show_transitions (CL);
      CL << "tp_elim "; tp_elim.show_transitions (CL);
    }

  // OK, now we have done the hard work of null state elimination. The rest is just book-keeping...

  // figure out which states are bifurcation targets
  set<int> bifurc_targets;
  for_const_contents (vector<int>, emit_and_bifurc_states, s)
    if (is_bifurc_type (state_type[*s]))
      {
	bifurc_targets.insert (bifurc[*s].l);
	bifurc_targets.insert (bifurc[*s].r);
      }

  // fill c2b_state, keeping emit states, bifurcation states & bifurcation targets
  c2b_state.clear();
  map<int,int> b2c_state;
  for (int sb = 0; sb < states(); ++sb)
    if (is_emit_type (state_type[sb]) || is_bifurc_type (state_type[sb]) || bifurc_targets.find(sb) != bifurc_targets.end())
      {
	b2c_state[sb] = c2b_state.size();
	c2b_state.push_back (sb);
      }

  // create Pair SCFG c
  Pair_CFG_scores c_cfg (c2b_state.size());
  c_cfg.name.clear();
  c_cfg.name << name << CFG_ELIM_SUFFIX;

  // copy over states from b to c
  for (int sc = 0; sc < c_cfg.states(); ++sc)
    {
      const int sb = c2b_state[sc];
      c_cfg.state_name[sc] = state_name[sb];

      const State_type s_type = state_type[sb];
      const bool s_is_null = !is_emit_type(s_type) && !is_bifurc_type(s_type);
      if (s_is_null && bifurc_targets.find(sb) == bifurc_targets.end())
	THROWEXPR ("Ulp, found a null state that's not a bifurcation target. Panicking");

      // copy everything over
      c_cfg.state_name[sc] = state_name[sb];
      c_cfg.state_type[sc] = state_type[sb];
      c_cfg.emit[sc] = emit[sb];

      // do bifurcations
      if (is_bifurc_type (s_type))
	c_cfg.bifurc[sc] = Bifurcation (b2c_state[bifurc[sb].l], b2c_state[bifurc[sb].r]);

      // do transitions
      for (int dc = 0; dc < c_cfg.states(); ++dc)
	{
	  const int db = c2b_state[dc];
	  const State_type d_type = state_type[db];
	  const bool d_is_null = !is_emit_type(d_type) && !is_bifurc_type(d_type);
	  if (!d_is_null)
	    c_cfg.transition (sc, dc) = Prob2Score (tp_elim.transition (sb, db));
	}
      if (!s_is_null)
	c_cfg.start[sc] = Prob2Score (tp_elim.transition (Grammar_state_enum::Start, sb));
      c_cfg.end[sc] = Prob2Score (tp_elim.transition (sb, Grammar_state_enum::End));

      // copy silly meta stuff too
      c_cfg.xlmeta_idx[sc] = xlmeta_idx[sb];
      c_cfg.xrmeta_idx[sc] = xrmeta_idx[sb];
      c_cfg.ylmeta_idx[sc] = ylmeta_idx[sb];
      c_cfg.ylmeta_idx[sc] = ylmeta_idx[sb];
    }

  // return
  return c_cfg;
}

// sample paths through null states (i.e. from a C-path, restore the B-path).
// returned state path is in b co-ordinates
vector<int> Pair_CFG_scores::sample_eliminated_states (const Transition_probs& tp_orig,  // b co-ords
						       const Transition_probs& tp_elim,  // b co-ords
						       const vector<int>& c2b_state,  // c->b
						       const Pair_CFG_scores& b_cfg,
						       const vector<int>& c_path)  // c co-ords
{
  // do c->b state index translation
  vector<int> b_path (c_path.size());
  for (int i = 0; i < (int) c_path.size(); ++i)
    b_path[i] = c_path[i] < 0 ? c_path[i] : c2b_state[c_path[i]];

  // get b's nonemit states
  const vector<int> nonemit = b_cfg.nonemit_states_unsorted();

  // break path into parse tree branches, resample each branch and re-concatenate
  const Pair_CFG_parse_tree b_tree = b_cfg.parse (b_path);
  vector<int> resampled_b_path;
  for_const_contents (vector<Pair_CFG_branch>, b_tree, b_branch)
    {
      const vector<int> resampled_b_branch_path = Transition_methods::sample_eliminated (tp_orig, tp_elim, nonemit, b_branch->path);
      resampled_b_path.insert (resampled_b_path.end(), resampled_b_branch_path.begin(), resampled_b_branch_path.end());
    }

  // return
  return resampled_b_path;
}

// sample null (no-sequence emitting) subtrees (i.e. from a B-path, restore the A-path).
// returned state path is in a co-ordinates
vector<int> Pair_CFG_scores::sample_eliminated_bifurcations (const vector<Prob>& p_empty,  // a co-ords
							     const vector<int>& b2a_state,  // b->a
							     const Pair_CFG_scores& a_cfg,
							     const Pair_CFG_scores& b_cfg,
							     const vector<int>& b_path)  // b co-ords
{
  // whenever we have a transition from a bifurcation-preceding null state to a state that is NOT the corresponding bifurcation,
  // this represents an empty subtree (or two empty subtrees, if it's a transition to the end state).
  //
  // sketch of algorithm:
  //  translate b->a;
  //  build @parse_tree (array of branches);
  //  resampled path = join (map (resample_branch($_.branch_path), @parse_tree))
  //
  // where...
  //
  //  resample_branch(branch_path) :-
  //    return map (resample_transition(path[$_-1], path[$_]), 1..@path-1);
  //
  //  resample_transition(i,j) :-
  //    if (i is NOT a bifurcation-preceding state):
  //     return i-->j-->...;
  //    else if (i is a proxy for the original bifurcation state, and j is the corresponding duplicated bifurcation state):
  //     return i-->...;
  //    else (i is a bifurcation-preceding state and j is NOT i's corresponding bifurcation state):
  //     let bif denote the true bifurcation.
  //     if j is the end state:
  //      return i-->(sample_null_tree(bif.left), sample_null_tree(bif.right))
  //     if j = bif.left AND ( j != bif.right OR coinflip=heads )
  //      return i-->(j-->..., sample_null_tree(bif.right));
  //     if j = bif.right:
  //      return i-->(sample_null_tree(bif.left), j-->...);

  vector<int> a_path;
  const Pair_CFG_parse_tree b_tree = b_cfg.parse (b_path);
  vector<list<int> > a_path_r (b_tree.size(), list<int>());  // stacks for R-emitted null subtrees; indexed by b_tree branch index
  for (int b_branch_idx = 0; b_branch_idx < (int) b_tree.size(); ++b_branch_idx)
    {
      const Pair_CFG_branch& b_branch = b_tree[b_branch_idx];
      const vector<int>& b_branch_path = b_branch.path;
      list<int>& a_path_r_branch = a_path_r[b_branch_idx];
      if (b_branch_path.size())
	{
	  const int b_path_start = b_branch_path.front();
	  a_path.push_back (b_path_start < 0 ? b_path_start : b2a_state[b_path_start]);
	  for (int pos = 1; pos < (int) b_branch_path.size(); ++pos)
	    {
	      const int ib = b_branch_path[pos - 1];
	      const int jb = b_branch_path[pos];
	      const int ia = ib < 0 ? ib : b2a_state[ib];
	      const int ja = jb < 0 ? jb : b2a_state[jb];
	      const bool i_is_bifurc_prec = ia < 0 ? false : (is_bifurc_type (a_cfg.state_type[ia]) && !is_bifurc_type (b_cfg.state_type[ib]));
	      if (!i_is_bifurc_prec)
		a_path.push_back (ja);  // normal transition
	      else // i_is_bifurc_prec
		{
		  const bool j_is_bifurc = ja == ia && jb != ib;
		  // if j_is_bifurc==true, then j is i's corresponding bifurcation state, i.e. a dummy state, and we don't add it to the a_path
		  if (!j_is_bifurc)
		    {
		      const Bifurcation& bif = a_cfg.bifurc[ia];
		      if (ja == End)
			{
			  // j is the End state, representing a double-null bifurcation; so add null L- and R-subtrees to a_path
			  const vector<int> l_null = sample_null_subtree (p_empty, a_cfg, bif.l);
			  const vector<int> r_null = sample_null_subtree (p_empty, a_cfg, bif.r);
			  a_path.insert (a_path.end(), l_null.begin(), l_null.end());
			  a_path.insert (a_path.end(), r_null.begin(), r_null.end());
			}
		      else if (ja == bif.l && (bif.l != bif.r || Rnd::decide(0.5)))
			{
			  // j is the left half of a bifurcation; so add j to a_path, and stack up a null R-subtree for when this branch & all its kids are finished
			  a_path.push_back (ja);
			  const vector<int> r_null = sample_null_subtree (p_empty, a_cfg, bif.r);
			  a_path_r_branch.insert (a_path_r_branch.begin(), r_null.begin(), r_null.end());
			}
		      else if (ja == bif.r)
			{
			  // j is the right half of a bifurcation; so add a null L-subtree to a_path, followed by j
			  const vector<int> l_null = sample_null_subtree (p_empty, a_cfg, bif.l);
			  a_path.insert (a_path.end(), l_null.begin(), l_null.end());
			  a_path.push_back (ja);
			}
		      else
			THROWEXPR ("Unknown transition: " << ib << "(" << ia << ") --> " << jb << "(" << ja << "); bifurcation is (" << bif.l << "," << bif.r << ")");
		    }
		}
	    }

	  // deal with stacked-up null R-subtrees: if we're childless they're our problem, otherwise they're inherited by our right child
	  if (b_branch.has_children())
	    a_path_r_branch.swap (a_path_r[b_branch.rchild]);
	  else
	    a_path.insert (a_path.end(), a_path_r_branch.begin(), a_path_r_branch.end());
	}
    }

  // return
  return a_path;
}
