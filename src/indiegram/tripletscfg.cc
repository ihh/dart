#include <stack>
#include "scfg/foldenv.h"
#include "indiegram/tripletscfg.h"
#include "util/strsaver.h"

// ancestral secondary structure
#define Stockholm_ancestral_secondary_structure_tag  "ancestral_SS"

// define these for speed
// used by SCFG_state_typing::_emit_size[]
#define _emit_N 0 // N as in null
#define _emit_0 1
#define _emit_1 SCFG_alphabet_size
#define _emit_2 (SCFG_alphabet_size * SCFG_alphabet_size)
#define _emit_3 (SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size)
#define _emit_4 (SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size)
#define _emit_5 (SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size)
#define _emit_6 (SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size * SCFG_alphabet_size)

sstring RNA_triplet_path::foldstring_ungapped_x() const
{
  sstring fold;
  for_const_contents (sstring, foldstring_x, c)
    if (*c != GAP_CHAR)
      fold.push_back (*c);
  return fold;
}

sstring RNA_triplet_path::foldstring_ungapped_y() const
{
  sstring fold;
  for_const_contents (sstring, foldstring_y, c)
    if (*c != GAP_CHAR)
      fold.push_back (*c);
  return fold;
}

sstring RNA_triplet_path::foldstring_ungapped_z() const
{
  sstring fold;
  for_const_contents (sstring, foldstring_z, c)
    if (*c != GAP_CHAR)
      fold.push_back (*c);
  return fold;
}

int RNA_triplet_path::foldstring_ungapped_x_len() const
{
  int len = 0;
  for_const_contents (sstring, foldstring_x, c)
    if (*c != GAP_CHAR)
      ++len;
  return len;
}

int RNA_triplet_path::foldstring_ungapped_y_len() const
{
  int len = 0;
  for_const_contents (sstring, foldstring_y, c)
    if (*c != GAP_CHAR)
      ++len;
  return len;
}

int RNA_triplet_path::foldstring_ungapped_z_len() const
{
  int len = 0;
  for_const_contents (sstring, foldstring_z, c)
    if (*c != GAP_CHAR)
      ++len;
  return len;
}

// fill this in (low priority?): see the .h file for an explanation of what it should do
Alignment_path RNA_triplet_path::alignment_path() const
{
  if ((foldstring_x.size() != foldstring_y.size()) || (foldstring_x.size() != foldstring_z.size()) || (foldstring_x.size() != foldstring_w.size()))
    THROWEXPR ("Gapped fold string have different lengths: X = '" << foldstring_x << "', Y = '" << foldstring_y << "', Z = '" << foldstring_z << "', W = '" << foldstring_w << "'");

  Alignment_path path;

  return path;

}

// constructor
Triplet_SCFG_alignment::Triplet_SCFG_alignment (const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z)
  : np_x (np_x), np_y (np_y), np_z (np_z)
{
  seqlen_x = np_x.size();
  seqlen_y = np_y.size();
  seqlen_z = np_z.size();
}

void Triplet_SCFG_alignment::show (ostream& o) const
{
  save_flags (o);
  left_align (o);

  if ((foldstring_x.size() != foldstring_y.size()) || (foldstring_x.size() != foldstring_z.size()) || (foldstring_x.size() != foldstring_w.size()))
    THROWEXPR ("Gapped fold string have different lengths: X = '" << foldstring_x << "', Y = '" << foldstring_y << "', Z = '" << foldstring_z << "', W = '" << foldstring_w << "'");

  sstring gappedw (foldstring_w.size());
  sstring gappedx (foldstring_x.size());
  sstring gappedy (foldstring_y.size());
  sstring gappedz (foldstring_z.size());
  int i, xpos, ypos, zpos;

  // read in the (gapped) aligned sequences
  // if foldstring_x[i] == GAP_CHAR then GAP_CHAR else np_x.seq[xpos++]
  for (i = 0, xpos = start_x; i < (int) foldstring_x.size(); ++i)
    gappedx[i] = foldstring_x[i] == GAP_CHAR ? GAP_CHAR : np_x.seq[xpos++];

  for (i = 0, ypos = start_y; i < (int) foldstring_y.size(); ++i)
    gappedy[i] = foldstring_y[i] == GAP_CHAR ? GAP_CHAR : np_y.seq[ypos++];

  for (i = 0, zpos = start_z; i < (int) foldstring_z.size(); ++i)
    gappedz[i] = foldstring_z[i] == GAP_CHAR ? GAP_CHAR : np_z.seq[zpos++];

  sstring ps_cons (foldstring_x.size()); // primary sequence conservation
  sstring ss_cons (foldstring_x.size()); // secondary structure conservation
  for (i = 0; i < (int) foldstring_x.size(); ++i)
    {
      ps_cons[i] = (gappedx[i] == gappedy[i] && gappedx[i] == gappedz[i]) ? gappedx[i] : '.';
      if (is_lchar(foldstring_x[i]) && is_lchar(foldstring_y[i]) && is_lchar(foldstring_z[i]))
	ss_cons[i] = FOLD_LCHAR;
      else if (is_rchar(foldstring_x[i]) && is_rchar(foldstring_y[i]) && is_rchar(foldstring_z[i]))
	ss_cons[i] = FOLD_RCHAR;
      else
	ss_cons[i] = NOFOLD_CHAR;
    }

  sstring xname, yname, zname, bits, xss, yss, zss, ps_cons_name, ss_cons_name, ancestral_ss_name;
  xname << np_x.name << '/' << start_x+1 << '-' << xpos;
  yname << np_y.name << '/' << start_y+1 << '-' << ypos;
  zname << np_z.name << '/' << start_z+1 << '-' << zpos;
  xss << Stockholm_sequence_column_annotation << ' ' << xname << ' ' << Stockholm_secondary_structure_tag;
  yss << Stockholm_sequence_column_annotation << ' ' << yname << ' ' << Stockholm_secondary_structure_tag;
  zss << Stockholm_sequence_column_annotation << ' ' << zname << ' ' << Stockholm_secondary_structure_tag;
  ps_cons_name << Stockholm_column_annotation << ' ' << StockholmConsensus(Stockholm_primary_sequence_tag);
  ss_cons_name << Stockholm_column_annotation << ' ' << StockholmConsensus(Stockholm_secondary_structure_tag);
  ancestral_ss_name << Stockholm_column_annotation << ' ' << Stockholm_ancestral_secondary_structure_tag;
  const int name_len = max (max( max( max(xss.size(),yss.size()),zss.size()), max (ss_cons_name.size(), bits.size())), ancestral_ss_name.size());

  o << "# " << Stockholm_format_identifier << " " << Stockholm_version_identifier << "\n";

  o.width (name_len);
  o << xss << " " << foldstring_x << "\n";
  o.width (name_len);
  o << xname << " " << gappedx << "\n";

  o.width (name_len);
  o << yss << " " << foldstring_y << "\n";
  o.width (name_len);
  o << yname << " " << gappedy << "\n";

  o.width (name_len);
  o << zss << " " << foldstring_z << "\n";
  o.width (name_len);
  o << zname << " " << gappedz << "\n";

  o.width (name_len);
  o << ancestral_ss_name << " " << foldstring_w << "\n";

  o.width (name_len);
  o << ps_cons_name << " " << ps_cons << "\n";

  o.width (name_len);
  o << ss_cons_name << " " << ss_cons << "\n";

  o.width (name_len);
  o << bits << " " << Score2Bits (score) << "\n";

  o << Stockholm_alignment_separator << "\n";

  restore_flags (o);

}

// constructor
Triplet_SCFG_branch::Triplet_SCFG_branch (int xl, int xr, int yl, int yr, int zl, int zr, int p)
  : xl (xl), xr (xr), yl (yl), yr (yr), zl (zl), zr (zr),
    parent (p), lchild (-1), rchild (-1)
{ }


// Recursive method to get the alignment corresponding to the subtree parse_tree.
// state_type is used to get the state type of each state
Triplet_SCFG_alignment Triplet_SCFG_branch::alignment (const vector<Triplet_SCFG_branch>& parse_tree, const vector<State_type>& state_type, const State_type_ancestral_map& state_type_ancestral_map, const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z) const
{

  sstring wlfold, xlfold, ylfold, zlfold, wrfold, xrfold, yrfold, zrfold;
  vector<State_type> parse;

  int xlpos = xl;
  int xrpos = xr;
  int ylpos = yl;
  int yrpos = yr;
  int zlpos = zl;
  int zrpos = zr;

  // go through all states in the state_path of the current branch
  // create an alignment as appropriate to the emission profiles of the emit states
  for (int i = 0; i < (int) (*this).path.size(); ++i)
    {
      const int s = (*this).path[i]; // state
      if (s < 0) continue; // ignore Start and End states
      const State_type t = state_type[s]; // state type

      State_type_ancestral_map::const_iterator stam_iter = state_type_ancestral_map.find (s);
      const bool is_emit_wl = stam_iter == state_type_ancestral_map.end() ? false : ((stam_iter->second & EmitWL) != 0);
      const bool is_emit_wr = stam_iter == state_type_ancestral_map.end() ? false : ((stam_iter->second & EmitWR) != 0);

      const bool is_emit_xl = ((t & EmitXL) != 0);
      const bool is_emit_xr = ((t & EmitXR) != 0);

      const bool is_emit_yl = ((t & EmitYL) != 0);
      const bool is_emit_yr = ((t & EmitYR) != 0);

      const bool is_emit_zl = ((t & EmitZL) != 0);
      const bool is_emit_zr = ((t & EmitZR) != 0);

      // emits to left or right?
      const bool is_emit_l = is_emit_wl || is_emit_xl || is_emit_yl || is_emit_zl;
      const bool is_emit_r = is_emit_wr || is_emit_xr || is_emit_yr || is_emit_zr;

      // symbols in sequences
      int xlc  = 0, xrc = 0, ylc = 0, yrc = 0, zlc = 0, zrc = 0;
      if (is_emit_xl) xlc = np_x.dsq[xlpos++];
      if (is_emit_xr) xrc = np_x.dsq[--xrpos];
      if (is_emit_yl) ylc = np_y.dsq[ylpos++];
      if (is_emit_yr) yrc = np_y.dsq[--yrpos];
      if (is_emit_zl) zlc = np_z.dsq[zlpos++];
      if (is_emit_zr) zrc = np_z.dsq[--zrpos];

      // emits basepairs?
      const bool is_paired_w = is_emit_wl && is_emit_wr;
      const bool is_paired_x = is_emit_xl && is_emit_xr;
      const bool is_paired_y = is_emit_yl && is_emit_yr;
      const bool is_paired_z = is_emit_zl && is_emit_zr;

      // emits Watson-Crick basepairs?
      const bool is_wcx = is_paired_x && (xlc == SCFG_alphabet.complement (xrc));
      const bool is_wcy = is_paired_y && (ylc == SCFG_alphabet.complement (yrc));
      const bool is_wcz = is_paired_z && (zlc == SCFG_alphabet.complement (zrc));

//      cerr << "i = " << i << "; s = " << s << "; t = " << t << "(" << state_type_string(t) << ")\n";
//      cerr << "is_emit_xl = " << is_emit_xl << "; is_paired_x = " << is_paired_x << "; is_wcx = " << is_wcx << "\n";
//      cerr << "is_emit_yl = " << is_emit_yl << "; is_paired_y = " << is_paired_y << "; is_wcy = " << is_wcy << "\n";
//      cerr << "is_emit_zl = " << is_emit_zl << "; is_paired_z = " << is_paired_z << "; is_wcz = " << is_wcz << "\n";
//      cerr << "xlpos = " << xlpos << "; xrpos = " << xrpos << "; ylpos = " << ylpos << "; yrpos = " << yrpos << "; zlpos = " << zlpos << "; zrpos = " << zrpos << "\n";
      

      // construct the local foldstring for this branch
      if (is_emit_l)
	{
	  // if is_emit_xl then { if is_paired_x then { if is_wcx then ... else ... } else ... } else ... }
	  wlfold << (is_emit_wl ? (is_paired_w ? (char) FOLD_LCHAR : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  xlfold << (is_emit_xl ? (is_paired_x ? (is_wcx ? (char) FOLD_LCHAR : (char) NONWC_LCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  ylfold << (is_emit_yl ? (is_paired_y ? (is_wcy ? (char) FOLD_LCHAR : (char) NONWC_LCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  zlfold << (is_emit_zl ? (is_paired_z ? (is_wcz ? (char) FOLD_LCHAR : (char) NONWC_LCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	}
      if (is_emit_r)
	{
	  wrfold << (is_emit_wr ? (is_paired_w ? (char) FOLD_RCHAR : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  xrfold << (is_emit_xr ? (is_paired_x ? (is_wcx ? (char) FOLD_RCHAR : (char) NONWC_RCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  yrfold << (is_emit_yr ? (is_paired_y ? (is_wcy ? (char) FOLD_RCHAR : (char) NONWC_RCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	  zrfold << (is_emit_zr ? (is_paired_z ? (is_wcz ? (char) FOLD_RCHAR : (char) NONWC_RCHAR) : (char) NOFOLD_CHAR) : (char) GAP_CHAR);
	}
      if (t != Null)
	{
	  parse.push_back (t);
	}
    }

  // if the current branch has children, recursively deal with them
  if ((*this).has_children())
    {
      // parse_tree[lchild] is the Triplet_SCFG_branch corresponding to the left child of the current branch
      const Triplet_SCFG_alignment lchild_alignment = parse_tree[lchild].alignment (parse_tree, state_type, state_type_ancestral_map, np_x, np_y, np_z);
      const Triplet_SCFG_alignment rchild_alignment = parse_tree[rchild].alignment (parse_tree, state_type, state_type_ancestral_map, np_x, np_y, np_z);
      xlfold << lchild_alignment.foldstring_x << rchild_alignment.foldstring_x; // by this point the (recursed) algorithm will have swapped the old foldstring_x for the new one
      ylfold << lchild_alignment.foldstring_y << rchild_alignment.foldstring_y; // it's easier to make the insertion on the left side (doing it on the right would require
      zlfold << lchild_alignment.foldstring_z << rchild_alignment.foldstring_z; // funny reversal of characters)
      wlfold << lchild_alignment.foldstring_w << rchild_alignment.foldstring_w;
      parse.insert (parse.end(), lchild_alignment.parse.begin(), lchild_alignment.parse.end());
      parse.insert (parse.end(), rchild_alignment.parse.begin(), rchild_alignment.parse.end());
    }
  // if not, then we're done!
  else
    {
      parse.push_back (Null); // Null in 'parse' indicates a leaf node
    }

  // So now we have used recursive calls along the tree to build new foldstrings, xlfold and xrfold.
  // All that remains is to put this information together.

  // Reverse the contents of xrfold.  ((*this).path is top->bottom along the branch, which corresponds to reading the right-emissions right->left
  // => need to reverse the contents of xrfold in order to get left->right reading of the right-emissions
  reverse (xrfold.begin(), xrfold.end()); // reverse(v.begin(), v.end()) reverses the contents of vector v
  reverse (yrfold.begin(), yrfold.end());
  reverse (zrfold.begin(), zrfold.end());
  reverse (wrfold.begin(), wrfold.end());

  // collate all the fold strings
  wlfold << wrfold;
  xlfold << xrfold;
  ylfold << yrfold;
  zlfold << zrfold;

  // create a new instance 'alignment' of Triplet_SCFG_alignment to hold the alignment we've made
  Triplet_SCFG_alignment alignment (np_x, np_y, np_z); 

  alignment.parse.swap (parse);

  // store the new foldstrings
  alignment.foldstring_w.swap (wlfold);
  alignment.foldstring_x.swap (xlfold);
  alignment.foldstring_y.swap (ylfold);
  alignment.foldstring_z.swap (zlfold);

  // store the starting positions of the local alignment
  alignment.start_x = xl;
  alignment.start_y = yl;
  alignment.start_z = zl;

  // ok, all done!
  return alignment;

}



// get the alignment corresponding to this parse tree
// relies on Triplet_SCFG_branch::alignment()
Triplet_SCFG_alignment Triplet_SCFG_parse_tree::alignment (const vector<State_type>& state_type, const SCFG_state_typing::State_type_ancestral_map& state_type_ancestral_map, const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z) const
{
  // if the parse tree is empty, return an alignment with Named_profile information only
  // Triplet_SCFG_parse_tree inherits from vector<Triplet_SCFG_branch>, so the size() method can be called to see if the parse tree is empty
  if (size() == 0)
    {
      Triplet_SCFG_alignment alignment (np_x, np_y, np_z);
      return alignment;
    }
  // use the method Triplet_SCFG_branch::alignment(parse_tree, state_type, npx, npy, npz) to create the alignment for this parse tree
  // Pass the parse tree as an argument and call alignment() on the root branch.
  return front().alignment (*this, state_type, state_type_ancestral_map, np_x, np_y, np_z);
}

// constructor
Triplet_alignment_with_tracebacks::Triplet_alignment_with_tracebacks (const Triplet_SCFG_alignment& alignment)
  : alignment (alignment)
{ }

int Triplet_SCFG_parse_tree::new_branch (int xl, int xr, int yl, int yr, int zl, int zr, int p)
{
  (*this).push_back (Triplet_SCFG_branch (xl, xr, yl, yr, zl, zr, p));  // push the new branch onto the Triplet_SCFG_parse_tree
  return size()-1; // return the index of the new branch
}
 
int Triplet_SCFG_parse_tree::new_lchild (int xl, int xr, int yl, int yr, int zl, int zr, int p)
{
  const int lc = new_branch (xl, xr, yl, yr, zl, zr, p);
  (*this)[p].lchild = lc;
  return lc; // return the index of the branch for the added lchild
}

int Triplet_SCFG_parse_tree::new_rchild (int xl, int xr, int yl, int yr, int zl, int zr, int p)
{
  const int rc = new_branch (xl, xr, yl, yr, zl, zr, p);
  (*this)[p].rchild = rc;
  return rc; // return the index of the branch for the added rchild
}


bool Triplet_SCFG_parse_tree::test_connections() const
{
  if (size()) // if parse tree not empty
    {
      const Triplet_SCFG_branch& root = front();
      if (root.parent != -1) return 0; // sanity check
      int n_children = 0;
      for (int i = 0; i < (int) size(); i++)
	{
	  const Triplet_SCFG_branch& branch = (*this)[i];
	  if (branch.xl > branch.xr || branch.yl > branch.yr || branch.zl > branch.zr) return 0; // sanity check for subsequence coordinates
	  if (branch.path.size() == 0) return 0; // can't have branches with no states
	  if (branch.lchild >= 0) // if has left child
	    {
	      if (branch.rchild < 0) return 0; // must have right child if has left child
	      if (branch.lchild >= (int) size() || branch.rchild >= (int) size()) return 0; // if the indices for the left or right child branches are outside the bounds of the parse_tree
	      if (branch.lchild == branch.rchild) return 0; // can't have the same index
	      const Triplet_SCFG_branch& lchild = (*this)[branch.lchild];
	      const Triplet_SCFG_branch& rchild = (*this)[branch.rchild];
	      if (lchild.parent != i || rchild.parent != i) return 0;
	      if (lchild.xr != rchild.xl) return 0; // child subsequences must join smoothly
	      if (lchild.yr != rchild.yl) return 0;
	      if (lchild.zr != rchild.zl) return 0;
	      n_children += 2;
	    }
	  else
	    {
	      if (branch.rchild >= 0) return 0; // can't have right child w/o left child
	    }
	}
      if (n_children != (int) size()-1) return 0;  // -1 because we didn't count the root
    }

  return 1;

}
 
bool Triplet_SCFG_parse_tree::test_global (const Digitized_biosequence& dsq_x, const Digitized_biosequence& dsq_y, const Digitized_biosequence& dsq_z) const
{
  if (size() == 0) return 0; // must be at least a root!
  const Triplet_SCFG_branch& root = front();
  if (root.path.size() == 0) return 0; // must have at least one state
  if (root.path[0] != Grammar_state_enum::Start) return 0; // first state must be Start
  if (root.xl != 0 || root.yl != 0 || root.zl != 0) return 0; // corresponding subsequences must be the entire sequence
  if (root.xr != (int) dsq_x.size() || root.yr != (int) dsq_y.size() || root.zr != (int) dsq_z.size()) return 0;

  return 1;
  
}

void Triplet_SCFG_parse_tree::show (ostream& o, const vector<State_type>* state_type /* = 0 */) const
{
  typedef pair<int,int> Branch_indent; // first = branch, second = index for indentation when displaying
  stack<Branch_indent> next_branch;
  next_branch.push (Branch_indent (0,0));
  while (!next_branch.empty()) // look at the whole subtree
    {
      const Branch_indent bi = next_branch.top(); 
      const Triplet_SCFG_branch& b = (*this)[bi.first]; // get the current branch 
      const int indent = bi.second;
      next_branch.pop();
      for (int i = 0; i < indent; ++i) // indent as appropriate
	o << ' ';
      o << "[" << b.xl << ".." << b.xr << "," << b.yl << ".." << b.yr << "," << b.zl << ".." << b.zr << "]"; // print subsequence indices corresponding to this branch of the parse tree
      for_const_contents (vector<int>, b.path, s) // (type, container, item) iterate over the state path for this branch
	{
	  o << " " << *s; // print state
	  if (*s == Grammar_state_enum::Start)
	    o << " (Start)";
	  else if (*s == Grammar_state_enum::End)
	    o << " (End)";
	  else if (state_type) // default argument is null
	    o << " (" << state_type_string ((*state_type)[*s]) << ")"; 
	}
      o << "\n";
      
      if (b.has_children())
	{
	  next_branch.push (Branch_indent (b.rchild, indent+1)); // indent children one space
	  next_branch.push (Branch_indent (b.lchild, indent+1));
	}
    }

}


// Map from emit state types (profiles) to the number of possible emissions.
// ( used in SCFG_state_typing::emit_size() )
// Not elegant, but fast.

int SCFG_state_typing::_emit_size[64] = 
  { _emit_N, _emit_1, _emit_1, _emit_2, // _emit_size[0] = _emit_N = 0 because 0 corresponds to states of type null
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_3, _emit_4, _emit_4, _emit_5, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_3, _emit_4, _emit_4, _emit_5, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_3, _emit_4, _emit_4, _emit_5, 
    _emit_3, _emit_4, _emit_4, _emit_5, 
    _emit_4, _emit_5, _emit_5, _emit_6 };

int SCFG_state_typing::_emit_xl_mul[64] =
  { _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_3, _emit_4, _emit_4, _emit_5 };

int SCFG_state_typing::_emit_xr_mul[64] =
  { _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0, 
    _emit_N, _emit_0, _emit_N, _emit_0 };

int SCFG_state_typing::_emit_yl_mul[64] =
  { _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_2, _emit_3, _emit_3, _emit_4 };

int SCFG_state_typing::_emit_yr_mul[64] =
  { _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1, 
    _emit_N, _emit_N, _emit_0, _emit_1 };

int SCFG_state_typing::_emit_zl_mul[64] =
  { _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_1, _emit_2, _emit_2, _emit_3 };

int SCFG_state_typing::_emit_zr_mul[64] =
  { _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2, 
    _emit_N, _emit_N, _emit_N, _emit_N, 
    _emit_0, _emit_1, _emit_1, _emit_2 };


// tricky: makes use of state typing to invert the emit hash used in Triplet_SCFG::path_emit_score() and testtripletscfg.cc.
sstring SCFG_state_typing::emit_hash_to_string (int s, int emit_hash) const
{
  int t = state_type[s];
  sstring emit_string;

  // operations are performed in the order required by the hash (xl yl zl zr yr xr) we use
  // characters are extracted from right to left
  if (t & EmitXR) { int value = (emit_hash % (SCFG_alphabet_size*emit_xr_mul(t))); int xrc = value / emit_xr_mul(t); emit_string << (char) toupper (alphabet().int2char(xrc)); emit_hash -= value; }
  else { emit_string << (char) Fold_char_enum::GAP_CHAR; }

  if (t & EmitYR) { int value = (emit_hash % (SCFG_alphabet_size*emit_yr_mul(t))); int yrc = value / emit_yr_mul(t); emit_string << (char) toupper (alphabet().int2char(yrc)); emit_hash -= value; }
  else { emit_string << (char) Fold_char_enum::GAP_CHAR; }

  if (t & EmitZR) { int value = (emit_hash % (SCFG_alphabet_size*emit_zr_mul(t))); int zrc = value / emit_zr_mul(t); emit_string << (char) toupper (alphabet().int2char(zrc)); emit_hash -= value; }
  else { emit_string << (char) Fold_char_enum::GAP_CHAR; }

  if (t & EmitZL) { int value = (emit_hash % (SCFG_alphabet_size*emit_zl_mul(t))); int zlc = value / emit_zl_mul(t); emit_string << (char) toupper (alphabet().int2char(zlc)); emit_hash -= value; }
  else { emit_string << (char) Fold_char_enum::GAP_CHAR; }

  if (t & EmitYL) { int value = (emit_hash % (SCFG_alphabet_size*emit_yl_mul(t))); int ylc = value / emit_yl_mul(t); emit_string << (char) toupper (alphabet().int2char(ylc)); emit_hash -= value; }
  else { emit_string << (char) Fold_char_enum::GAP_CHAR; }

  if (t & EmitXL) { int value = (emit_hash % (SCFG_alphabet_size*emit_xl_mul(t))); int xlc = value / emit_xl_mul(t); emit_string << (char) toupper (alphabet().int2char(xlc)); emit_hash -= value; }
  else { emit_string << (char) Fold_char_enum::GAP_CHAR; }

  reverse (emit_string.begin(), emit_string.end());
  return emit_string;
}

// assignment operator
SCFG_state_typing& SCFG_state_typing::operator= (const SCFG_state_typing& base)
{
  assign_state_typing (base);
  return *this;
}

// state typing assignment
void SCFG_state_typing::assign_state_typing (const SCFG_state_typing& base)
{
  state_type = base.state_type;
  bifurc = base.bifurc;
  state_type_ancestral_map = base.state_type_ancestral_map;
}

vector<vector<SCFG_state_type_enum::Bifurcation_left_parent> > SCFG_state_typing::left_parent() const
{
  vector<vector<Bifurcation_left_parent> > lp (num_states());
  for (int s = 0; s < num_states(); ++s)
    if (is_bifurc_type (state_type[s]))
      lp[bifurc[s].r].push_back (Bifurcation_left_parent (bifurc[s].l, s));
  return lp;
}


vector<vector<SCFG_state_type_enum::Bifurcation_right_parent> > SCFG_state_typing::right_parent() const
{
  vector<vector<Bifurcation_right_parent> > rp (num_states());
  for (int s = 0; s < num_states(); ++s)
    if (is_bifurc_type (state_type[s]))
      rp[bifurc[s].l].push_back (Bifurcation_right_parent (bifurc[s].r, s));
  return rp;
}

const char* SCFG_state_type_enum::state_type_string (State_type t)
{
  switch (t)
    {
    case Null: return "Null"; break;
      // emit states
    case EmitXR: return "EmitXR"; break;
    case EmitYR: return "EmitYR"; break;
    case EmitXRYR: return "EmitXRYR"; break;
    case EmitZR: return "EmitZR"; break;
    case EmitXRZR: return "EmitXRZR"; break;
    case EmitYRZR: return "EmitYRZR"; break;
    case EmitXRYRZR: return "EmitXRYRZR"; break;
    case EmitZL: return "EmitZL"; break;
    case EmitXRZL: return "EmitXRZL"; break;
    case EmitYRZL: return "EmitYRZL"; break;
    case EmitXRYRZL: return "EmitXRYRZL"; break;
    case EmitZLR: return "EmitZLR"; break;
    case EmitXRZLR: return "EmitXRZLR"; break;
    case EmitYRZLR: return "EmitYRZLR"; break;
    case EmitXRYRZLR: return "EmitXRYRZLR"; break;
    case EmitYL: return "EmitYL"; break;
    case EmitXRYL: return "EmitXRYL"; break;
    case EmitYLR: return "EmitYLR"; break;
    case EmitXRYLR: return "EmitXRYLR"; break;
    case EmitYLZR: return "EmitYLZR"; break;
    case EmitXRYLZR: return "EmitXRYLZR"; break;
    case EmitYLRZR: return "EmitYLRZR"; break;
    case EmitXRYLRZR: return "EmitXRYLRZR"; break;
    case EmitYLZL: return "EmitYLZL"; break;
    case EmitXRYLZL: return "EmitXRYLZL"; break;
    case EmitYLRZL: return "EmitYLRZL"; break;
    case EmitXRYLRZL: return "EmitXRYLRZL"; break;
    case EmitYLZLR: return "EmitYLZLR"; break;
    case EmitXRYLZLR: return "EmitXRYLZLR"; break;
    case EmitYLRZLR: return "EmitYLRZLR"; break;
    case EmitXRYLRZLR: return "EmitXRYLRZLR"; break;
    case EmitXL: return "EmitXL"; break;
    case EmitXLR: return "EmitXLR"; break;
    case EmitXLYR: return "EmitXLYR"; break;
    case EmitXLRYR: return "EmitXLRYR"; break;
    case EmitXLZR: return "EmitXLZR"; break;
    case EmitXLRZR: return "EmitXLRZR"; break;
    case EmitXLYRZR: return "EmitXLYRZR"; break;
    case EmitXLRYRZR: return "EmitXLRYRZR"; break;
    case EmitXLZL: return "EmitXLZL"; break;
    case EmitXLRZL: return "EmitXLRZL"; break;
    case EmitXLYRZL: return "EmitXLYRZL"; break;
    case EmitXLRYRZL: return "EmitXLRYRZL"; break;
    case EmitXLZLR: return "EmitXLZLR"; break;
    case EmitXLRZLR: return "EmitXLRZLR"; break;
    case EmitXLYRZLR: return "EmitXLYRZLR"; break;
    case EmitXLRYRZLR: return "EmitXLRYRZLR"; break;
    case EmitXLYL: return "EmitXLYL"; break;
    case EmitXLRYL: return "EmitXLRYL"; break;
    case EmitXLYLR: return "EmitXLYLR"; break;
    case EmitXLRYLR: return "EmitXLRYLR"; break;
    case EmitXLYLZR: return "EmitXLYLZR"; break;
    case EmitXLRYLZR: return "EmitXLRYLZR"; break;
    case EmitXLYLRZR: return "EmitXLYLRZR"; break;
    case EmitXLRYLRZR: return "EmitXLRYLRZR"; break;
    case EmitXLYLZL: return "EmitXLYLZL"; break;
    case EmitXLRYLZL: return "EmitXLRYLZL"; break;
    case EmitXLYLRZL: return "EmitXLYLRZL"; break;
    case EmitXLRYLRZL: return "EmitXLRYLRZL"; break;
    case EmitXLYLZLR: return "EmitXLYLZLR"; break;
    case EmitXLRYLZLR: return "EmitXLRYLZLR"; break;
    case EmitXLYLRZLR: return "EmitXLYLRZLR"; break;
    case EmitXLRYLRZLR: return "EmitXLRYLRZLR"; break;
      // end emit states
    case Bifurc: return "Bifurc"; break;
    case Undefined: return "Undefined"; break;
    default: return "[unknown]"; break;
    }
  return "UnreachableText";
}

// to do: get rid of zero_outgoing_bifurcation_transitions()
// call Transition_methods::selected_outgoing_states() and then 
// add in bifurcation states by hand
vector<vector<int> > Triplet_SCFG::selected_outgoing_states (const vector<int>& selection) const
{
  // hack to avoid returning outgoing bifurcation states:
  //  set all transitions to bifurcation states to -InfinityScore
  //  so that Transition_methods::selected_outgoing_states () doesn't return them
  ((Triplet_SCFG*) this)->zero_outgoing_bifurcation_transitions();  // cast away const.. hacky
  return Transition_methods::selected_outgoing_states ((*this).transition_scores, selection);
}

// to do: fix this also
vector<vector<int> > Triplet_SCFG::selected_incoming_states (const vector<int>& selection) const
{
  ((Triplet_SCFG*) this) -> zero_outgoing_bifurcation_transitions();  // cast away const.. hacky
  return Transition_methods::selected_incoming_states ((*this).transition_scores, selection);
}


void Triplet_SCFG::zero_outgoing_bifurcation_transitions()
{
  for (int s = 0; s < num_states(); ++s)
    if (is_bifurc_type (state_type[s]))
      for (int d = 0; d < num_states(); ++d)
	transition_scores.transition (s, d) = -InfinityScore;
}


void Triplet_SCFG::add_fake_bifurcation_transitions()
{
  zero_outgoing_bifurcation_transitions();
  for (int s = 0; s < num_states(); ++s)
    if (is_bifurc_type (state_type[s]))
      {
	transition_scores.transition (s, bifurc[s].l) = 0;
	transition_scores.transition (s, bifurc[s].r) = 0;
      }
}


// minimal constructor
Triplet_SCFG::Triplet_SCFG()
  : SCFG_state_typing(),
    transition_scores (Transition_scores()), // NB: No constructor Concrete_transition_scores() constructor, so use Transition_scores()
    emit()
{ }


// constructor
Triplet_SCFG::Triplet_SCFG (int num_states)
  : SCFG_state_typing (num_states),
    transition_scores (Concrete_transition_scores (num_states, -InfinityScore)),
    emit (num_states)
{ }

// copy constructor
Triplet_SCFG::Triplet_SCFG (const Triplet_SCFG& base)
  : SCFG_state_typing (base),
    transition_scores (base.transition_scores),
    emit (base.emit)
{ }

// assignment operator
Triplet_SCFG& Triplet_SCFG::operator= (const Triplet_SCFG& scfg)
{
  assign_state_typing (scfg);
  transition_scores.assign_transition_matrix (scfg.transition_scores);
  emit = scfg.emit;
  return *this;
}

// code for the switch over Emit states created with this horrible one-liner on the state type list  Pair_CFG_state_typing:
// perl -a -ne'foreach my $f (@F) { if ($f =~ /^Emit[a-zA-Z]*$/) {print "case Pair_CFG_scores::$f:\ninit_emit (s, $f);\nfor (int xl = 0; xl < 4; ++xl)\nfor(int xr = 0; xr < 4; ++xr)\nfor(int yl = 0; yl < 4; ++yl)\nfor(int yr = 0; yr < 4; ++yr)\n{\n int triplet_emit_idx = emit_idx ($f, xl, xr, yl, yr, 0, 0);\nint pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::$f)\n + xr * paircfg.emit_xr_mul (Pair_CFG_scores::$f)\n + yl * paircfg.emit_yl_mul (Pair_CFG_scores::$f)\n + yr * paircfg.emit_yr_mul (Pair_CFG_scores::$f); \nemit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];\n}\nbreak;\n";} };'
Triplet_SCFG::Triplet_SCFG (const Pair_CFG_scores& paircfg)
  : SCFG_state_typing (paircfg.states()),
    transition_scores (paircfg.states()),
    emit (paircfg.states())
{
  if (&alphabet() != &paircfg.alphabet()) THROWEXPR ("Alphabet mismatch");
  this->transition_scores.start_to_end() = paircfg.start_to_end();
  for (int s = 0; s < num_states(); ++s)
    {
      // handle transition from Start and to End
      this->transition_scores.start[s] = paircfg.start[s];
      this->transition_scores.end[s] = paircfg.end[s];
      // handle all other transitions
      for (int d = 0; d < num_states(); ++d)
	this->transition_scores.transition (s, d) = paircfg.transition (s, d);
      // initialize States as appropriate
      switch (paircfg.state_type[s])
	{

	case Pair_CFG_scores::EmitXL:
	  init_emit (s, EmitXL);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXL, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXL)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXL)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXL)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXL); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXR:
	  init_emit (s, EmitXR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLR:
	  init_emit (s, EmitXLR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitYL:
	  init_emit (s, EmitYL);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitYL, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitYL)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitYL)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitYL)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitYL); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitYR:
	  init_emit (s, EmitYR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitYR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitYR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitYR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitYR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitYR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitYLR:
	  init_emit (s, EmitYLR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitYLR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitYLR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitYLR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitYLR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitYLR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLYL:
	  init_emit (s, EmitXLYL);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLYL, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLYL)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLYL)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLYL)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLYL); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLYR:
	  init_emit (s, EmitXLYR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLYR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLYR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLYR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLYR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLYR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLYLR:
	  init_emit (s, EmitXLYLR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLYLR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLYLR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLYLR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLYLR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLYLR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXRYL:
	  init_emit (s, EmitXRYL);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXRYL, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXRYL)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXRYL)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXRYL)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXRYL); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXRYR:
	  init_emit (s, EmitXRYR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXRYR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXRYR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXRYR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXRYR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXRYR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXRYLR:
	  init_emit (s, EmitXRYLR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXRYLR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXRYLR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXRYLR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXRYLR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXRYLR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLRYL:
	  init_emit (s, EmitXLRYL);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLRYL, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLRYL)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLRYL)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLRYL)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLRYL); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLRYR:
	  init_emit (s, EmitXLRYR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLRYR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLRYR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLRYR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLRYR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLRYR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;
	case Pair_CFG_scores::EmitXLRYLR:
	  init_emit (s, EmitXLRYLR);
	  for (int xl = 0; xl < 4; ++xl)
	    for(int xr = 0; xr < 4; ++xr)
	      for(int yl = 0; yl < 4; ++yl)
		for(int yr = 0; yr < 4; ++yr)
		  {
		    int triplet_emit_idx = emit_idx (EmitXLRYLR, xl, xr, yl, yr, 0, 0);
		    int pair_emit_idx = xl * paircfg.emit_xl_mul (Pair_CFG_scores::EmitXLRYLR)
		      + xr * paircfg.emit_xr_mul (Pair_CFG_scores::EmitXLRYLR)
		      + yl * paircfg.emit_yl_mul (Pair_CFG_scores::EmitXLRYLR)
		      + yr * paircfg.emit_yr_mul (Pair_CFG_scores::EmitXLRYLR); 
		    emit[s][triplet_emit_idx] = paircfg.emit[s][pair_emit_idx];
		  }
	  break;

	case Pair_CFG_scores::Null:
	  init_emit (s, Null);
	  break;
	case Pair_CFG_scores::Bifurc:
	  init_bifurc (s, paircfg.bifurc[s].l, paircfg.bifurc[s].r);
	  break;
	case Pair_CFG_scores::BifurcRevY:
	  THROWEXPR ("BifurcRevY not implemented for Triplet_SCFGs.  Offending state '" << s << "' with state_type[" << s << "] = '" << paircfg.state_type[s] << "'\n");
	default:
 	  THROWEXPR ("State type unknown for state '" << s << "' with state_type[" << s << "] = '" << paircfg.state_type[s] << "'\n");
	}
    }
}

// initialize an emit state
void Triplet_SCFG::init_emit (int state, State_type type)
{
  state_type[state] = type;
  emit[state] = vector<Score> (emit_size (state));
}

// initialize a bifurcation state
void Triplet_SCFG::init_bifurc (int state, int l, int r)
{
  state_type[state] = Bifurc;
  bifurc[state] = SCFG_state_type_enum::Bifurcation (l, r);
}


void Triplet_SCFG::indent_xr (Triplet_SCFG_parse_tree& parse_tree, int node, int offset)
{
  if (node < 0) return; // sanity check
  parse_tree[node].xr += offset;
  indent_xr (parse_tree, parse_tree[node].lchild, offset);
  indent_xr (parse_tree, parse_tree[node].rchild, offset);
}

void Triplet_SCFG::indent_yr (Triplet_SCFG_parse_tree& parse_tree, int node, int offset)
{
  if (node < 0) return; // sanity check
  parse_tree[node].yr += offset;
  indent_yr (parse_tree, parse_tree[node].lchild, offset);
  indent_yr (parse_tree, parse_tree[node].rchild, offset);
}

void Triplet_SCFG::indent_zr (Triplet_SCFG_parse_tree& parse_tree, int node, int offset)
{
  if (node < 0) return; // sanity check
  parse_tree[node].zr += offset;
  indent_zr (parse_tree, parse_tree[node].lchild, offset);
  indent_zr (parse_tree, parse_tree[node].rchild, offset);
}


// show method
void Triplet_SCFG::show (ostream& o) const
{
  int old_prec = o.precision(3);
  transition_scores.save_flags (o); // strsaver method (inherited by class Transition_matrix)
  transition_scores.right_align (o);

  transition_scores.show_transitions (o); // Transition_matrix method

  // to do: this doesn't actually show the scores, but rather an integer indexing the emission
  // what's the point of keeping it?  maybe get rid of it
  o << "Emission profile " << transition_scores.element_descriptor() << ":\n"; // 

  const int w = transition_scores.element_width() + 1;
  int max_sz = 0;  // find the maximum number of emissions from a state
  for (int s = 0; s < transition_scores.tm_states(); ++s) { max_sz = max (max_sz, emit_size (s)); }
  for (; max_sz > 1; max_sz = max_sz / alphabet_size())
    {
      o << "                     ";
      for (int s = 0; s < max_sz; s++)
	{
	  sstring text;
	  for (int i = 1; i < max_sz; i *= alphabet_size())
	    text << (char) toupper (alphabet().int2char((s/i) % alphabet_size()));
	  text << "=" << s;
	  o.width (w+1);
	  o << text;
	}
      o << "\n";
    }

  // now show the emissions for all states as well as the associated emit scores
  for (int s = 0; s < transition_scores.tm_states(); s++)
    {
      sstring text;
      text << "State " << s << "[" << state_type_string (state_type[s]) << "]"; // State [state_type]
      transition_scores.left_align (o); // Transition_matrix inherits from Stream_saver
      o.width (20);
      o << text << " (";
      transition_scores.right_align (o);

      // for all possible emissions x of state s
      for (int x = 0; x < (int) emit[s].size(); x++)
	{
	  o.width (w);
	  o << emit_hash_to_string (s, x) << "=";
	  show_element (emit[s][x], o);
	  if (x < (int) emit[s].size() - 1) // for prettification
	    o << " ";
	}
      o << ")";

      //      for_const_contents (vector<int>, 

      // show bifurcations
      if (is_bifurc_type (state_type[s]))      // if (state_type[s] == Bifurc)
	o << " (bifurcates to " << bifurc[s].l << ", " << bifurc[s].r << ")";
      o << "\n";
    }

  transition_scores.restore_flags (o);
  o.precision (old_prec);

}

// state_path corresponds to a depth-first, left->right traversal of the parse tree
Triplet_SCFG_parse_tree Triplet_SCFG::parse (const vector<int>& state_path, int xstart, int ystart, int zstart) const
{
  // xl, xr keep track of the number of residues emitted to X to the left, right along the current branch
  // Remember that Triplet_SCFG_branch has members xl, xr also which store the coordinates of the subsequence of X
  // it represents.
  int xl = 0;
  int xr = 0;
  int yl = 0;
  int yr = 0;
  int zl = 0;
  int zr = 0;
  
  // xlen[node] is the (integer) number of symbols emitted to X by the current branch and all of its children (the subtree rooted at X)
  // it's used to construct the appropriate right subsequence indices
  vector<int> xlen (1, (int) 0);
  vector<int> ylen (1, (int) 0);
  vector<int> zlen (1, (int) 0);

  // initialize parse tree (create root branch with 3 subsequences of length 0)
  Triplet_SCFG_parse_tree parse_tree (0, 0, 0);

  // Build the parse tree: loop over all states in the state path.
  for (int i = 0; i < (int) state_path.size(); ++i)
    {
      const int s = state_path[i]; // get current state s
      //cerr << "i = " << i << "; state s = " << s << "; type t = " << t << "\n";
      // Excepting the root branch, branches begin at bifurcations and end at End states.
      // For each branch, we need to keep track of: the indices of the corresponding X, Y, Z subsequences,
      // the parse tree connections (meaning the parent), and the state path.
      // We therefore keep track of:
      // - length of emissions to the left and right along the current branch (xl, xr)
      // - total length (left and right concatenated) of emissions along the subtree rooted at the current branch (xlen, ylen):
      //    this is necessary for updating the right subsequence indices (we don't a priori know the total sequence lenths represented by the parse tree,
      //    so we have to figure them out as we go
      // - the parent node of the current branch
      // - the state path.

      // get the branch which we're working on right now
      Triplet_SCFG_branch& current_branch = parse_tree.back();

      // add the state which we're looking at now to the state path of the current branch
      current_branch.path.push_back (s);

      // sanity checks: only the first state in state_path can be the Start state
      if (i == 0) // had better be the Start state!
	{
	  if (s != Grammar_state_enum::Start) THROWEXPR ("No Start state");  // (note that Start, End are handled as special cases; they're not state types)
	  continue; // nothing more to do with the start state
	}
      if (s == Grammar_state_enum::Start) THROWEXPR ("Unexpected Start state");



      // Now build the tree by constructing branches and connecting them as appropriate.

      // if we're in the middle of creating a branch
      if (s != Grammar_state_enum::End)
	{
	  const State_type t = state_type[s];

	  // sanity check: if this is the last state in the state path, then something's wrong: it should be an End state
	  if (i == (int) state_path.size() - 1) THROWEXPR ("Unterminated branch");

	  // if a bifurcation state
	  if (is_bifurc_type (t)) // if (t == Bifurc) // this test should be equivalent to if (is_bifurc_type (t)); Triplet_SCFG inherits member state_type from SCFG_state_typing
	    {
	      // hit a bifurcation, so create the corresponding branch p and stick it at the end of the parse tree built so far
	      const int p = parse_tree.size() - 1;
	      //cerr << "bifurcation terminates branch " << p << ": s = " << s << "\n";
	      // we're doing a depth-first, left->right traversal of the parse tree, so the next state will lie on the left child branch of the parent branch p:
	      // create a left child of the bifurcation branch
	      // and update the branch subsequence indices as appropriate given the emissions so far along the current branch.
	      parse_tree.new_lchild (current_branch.xl + xl, current_branch.xr + xr, current_branch.yl + yl, current_branch.yr + yr, current_branch.zl + zl, current_branch.zr + zr, p);
	      xlen.push_back (0); // create emission entries for the newly created branch and set to zero
	      ylen.push_back (0);
	      zlen.push_back (0);
	      // we're on to a new branch, so reset the emission bookkeeping indices
	      xl = xr = yl = yr = zl = zr = 0;
	    }

	  // if an emit or null state
	  else
	    {
	      // if state s left-emits to X, increment the emission bookkeeping indices as appropriate
	      if (t & EmitXL) { ++xl; ++xlen.back(); }
	      if (t & EmitXR) { ++xr; ++xlen.back(); }
	      if (t & EmitYL) { ++yl; ++ylen.back(); }
	      if (t & EmitYR) { ++yr; ++ylen.back(); }
	      if (t & EmitZL) { ++zl; ++zlen.back(); }
	      if (t & EmitZR) { ++zr; ++zlen.back(); }
	    }
	}
      // if we've just finished a branch, perform the necessary bookkeeping
      else
	{
	  // recurse up the tree, looking for the last branch without a right child
	  // also increment the subtree emission indices xlen, etc. as we go and update the branch subsequence indices as appropriate
	  int c = parse_tree.size() - 1; // child index; start at the last branch we worked on
	  int p = -1;
	  //cerr << "finished a branch: c = " << c << "; p = " << p << "; parse_tree[p].rchild = " << parse_tree[p].rchild << "\n";
	  while (c != 0) // note that this recursion /does/ properly update the root branch (c != 0 makes sure that it doesn't try to update the parent of the root branch)
	    {
	      p = parse_tree[c].parent;
	      // increment subtree emission indices
	      xlen[p] += xlen[c];
	      ylen[p] += ylen[c];
	      zlen[p] += zlen[c];

	      //cerr << "c = " << c << "; p = " << p << "; parse_tree[p].rchild = " << parse_tree[p].rchild << "\n";
	      // if we encounter a parent branch without a right child, break and deal with it
	      if (parse_tree[p].rchild < 0) break;

	      // update the branch subsequence indices as we recurse up the tree
	      const int lc = parse_tree[p].lchild;
	      indent_xr (parse_tree, lc, xlen[c]);
	      indent_yr (parse_tree, lc, ylen[c]);
	      indent_zr (parse_tree, lc, zlen[c]);
	      c = p; // go up the tree one level
	    }

	  // we've located the deepest branch without a right child
	  // again, this /does/ properly catch the root branch case since it tests (c > 0), not (p > 0)
	  if (c > 0)
	    {
	      // sanity check: we found a branch without a right child, but we've run out of states in the state path to create the right child branch with!
	      if (i == (int) state_path.size() - 1) THROWEXPR ("Nonempty stack");
	      // create right child
	      const Triplet_SCFG_branch& left = parse_tree[c];
	      parse_tree.new_rchild (left.xl + xlen[c], left.xr, left.yl + ylen[c], left.yr, left.zl + zlen[c], left.zr, p);
	      xlen.push_back (0); // create emission indices for the newly created branch and set to 0
	      ylen.push_back (0);
	      zlen.push_back (0);
	    }
	  // sanity check
	  else // if (c == 0)
	    // sanity check: we've reached the root branch without encountering a branch without a right child,
	    // so we've finished constructing the parse tree and there shouldn't remain any states in the state path!
	    if (i < (int) state_path.size() - 1) { cerr << "c  = " << c << "; p = " << p << "; i = " << i << "; state_path.size() = " << state_path.size() << "\n"; THROWEXPR("Empty stack in Triplet_SCFG::parse()."); }
	
	  // finished a branch, so reset emission indices
	  xl = xr = yl = yr = zl = zr = 0;
	}

    }
  // we've finished building the parse tree from the state path

  // now update the branch subsequence coordinates per the xstart, ystart, zstart passed as parameters
  const int total_xlen = xlen[0];
  const int total_ylen = ylen[0];
  const int total_zlen = zlen[0];

  for_contents (Triplet_SCFG_parse_tree, parse_tree, branch)
    {
      branch->xl += xstart;
      branch->yl += ystart;
      branch->zl += zstart;
      branch->xr = xstart + total_xlen - branch->xr;
      branch->yr = ystart + total_ylen - branch->yr;
      branch->zr = zstart + total_zlen - branch->zr;
    }
  

  return parse_tree;
}

Triplet_SCFG_parse_tree Triplet_SCFG::parse (const vector<int>& state_path) const
{
  return parse (state_path, 0, 0, 0);
}

Score Triplet_SCFG::path_transition_score (const Triplet_SCFG_parse_tree& parse_tree) const
{
  Score score = 0;
  // Transition_methods::path_transition_score() is:
  // static Score path_transition_score (const Transition_scores& ts, const vector<int>& state_path);
  // where Transition_scores is a typedef for Transition_matrix<Score>
  // b->path is the state path along the branch b of parse_tree.
  for_const_contents (Triplet_SCFG_parse_tree, parse_tree, b)
    ScorePMulAcc (score, Transition_methods::path_transition_score ((*this).transition_scores, b->path));
  return score;
}

Score Triplet_SCFG::path_emit_score (const Triplet_SCFG_parse_tree& parse_tree, const Digitized_biosequence& dsq_x, const Digitized_biosequence& dsq_y, const Digitized_biosequence& dsq_z) const
{
  Score score = 0;
  for_const_contents (Triplet_SCFG_parse_tree, parse_tree, b)
    {
      // get sequence coordinates
      int xl = b->xl;
      int xr = b->xr;
      int yl = b->yl;
      int yr = b->yr;
      int zl = b->zl;
      int zr = b->zr;

      for_const_contents (vector<int>, b->path, s)
	if (*s >= 0) // for states != Start, End (no associated emit score with these guys)
	  {
	    const State_type t = state_type[*s];
	    if (is_emit_type (t))
	      {
		// get emit_idx and then increment the indices as appropriate (if an emission)
		int emit_idx = (*this).emit_idx (t, dsq_x[xl], dsq_x[xr-1], dsq_y[yl], dsq_y[yr-1], dsq_z[zl], dsq_z[zr-1]);
		if (t & EmitXL) ++xl;
		if (t & EmitXR) --xr;
		if (t & EmitYL) ++yl;
		if (t & EmitYR) --yr;
		if (t & EmitZL) ++zl;
		if (t & EmitZR) --zr;
		const Score emit_score = emit[*s][emit_idx];
		// debugging
		//cout << "state " << *s << " (type " << t << "): " << emit_hash_to_string (*s, emit_idx) << " = " << emit_score << "\n";
		//
		ScorePMulAcc (score, emit_score);
	      }
	  }
    }

  return score;
}

Score Triplet_SCFG::path_score (const Triplet_SCFG_parse_tree& parse_tree, const Digitized_biosequence& dsq_x, const Digitized_biosequence& dsq_y, const Digitized_biosequence& dsq_z) const
{
  return ScorePMul (path_transition_score (parse_tree), path_emit_score (parse_tree, dsq_x, dsq_y, dsq_z));
}
