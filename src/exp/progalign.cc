#include "handel/progalign.h"
#include "telegraph/tgvar.h"
#include "hmm/pairhmmem.h"
#include "hmm/pairenv.h"
#include "hmm/postpairhmm.h"
#include "util/vector_output.h"

// #define's

#define OPTACC_CONSISTENCY_ITERATIONS 2
#define OPTACC_REFINEMENT_ITERATIONS 100

// Progressive_aligner

Tree_alignment Progressive_aligner::make_progressive_alignment (const Sequence_database_index& seqdb_index)
{
  // initialise
  CTAG(5,PROGALIGN) << "Initialising progressive aligner\n";
  initialise (seqdb_index);

  // make guide tree
  CTAG(5,PROGALIGN) << "Making guide tree\n";
  make_guide_tree();
  if (CTAGGING(5,PROGALIGN))
    {
      CL << "Guide tree:\n";
      tree.write (CL);
    }

  // progressive phase
  CTAG(5,PROGALIGN) << "Doing progressive alignment\n";
  for_rooted_nodes_post (tree, branch_iter)
    {
      const Node parent = (*branch_iter).first;
      const Node node = (*branch_iter).second;
      const vector<Node> child (tree.children_begin(node,parent), tree.children_end(node,parent));
      if (child.size() == 0)
	continue;
      else if (child.size() == 2)
	{
	  // get x & y
	  const Node x = child[0];
	  const Node y = child[1];

	  // get tx & ty
	  const double tx = tree.branch_length (node, x);
	  const double ty = tree.branch_length (node, y);

	  // get names
	  const sstring aname = tree.node_specifier (node);
	  const sstring xname = tree.node_specifier (x);
	  const sstring yname = tree.node_specifier (y);
	  
	  // print log message
	  CTAG(4,PROGALIGN) << "Aligning '" << xname << "' (distance " << tx << ") and '" << yname << "' (distance " << ty << ") to get '" << aname << "'\n";

	  // do it
	  const Alignment_path axy_path = align_axy (node, x, y, tx, ty);
	  decomp[Row_pair (node, x)] = Pairwise_path (axy_path, 0, 1, TRUE);
	  decomp[Row_pair (node, y)] = Pairwise_path (axy_path, 0, 2, TRUE);
	  
	  // print log message
	  vector<sstring> row_name (3);
	  row_name[0] = tree.node_specifier (node);
	  row_name[1] = tree.node_specifier (x);
	  row_name[2] = tree.node_specifier (y);
	  if (CTAGGING(4,PROGALIGN))
	    {
	      CL << "Alignment of node '" << row_name[0] << "' and descendants:\n";
	      axy_path.show (CL, row_name);
	    }
	}
      else
	THROWEXPR ("Not a binary tree");
    }

  // print log message
  if (CTAGGING(3,PROGALIGN))
    {
      Alignment_path path;
      path.compose (decomp, TRUE);
      vector<sstring> dummy_row_name (path.rows());
      CTAG(3,PROGALIGN) << "Alignment path after progressive alignment:\n";
      path.show (CL, dummy_row_name);
    }

  // refinement phase
  CTAG(5,PROGALIGN) << "Doing iterative refinement\n";
  refine();

  // print log message
  CTAG(5,PROGALIGN) << "Completed alignment\n";
  // make Tree_alignment
  Tree_alignment tree_align;
  tree_align.set_tree (tree);
  tree_align.make_empty_alignment();
  for (int i = 0; i < seqdb_index.size(); ++i)
    tree_align.align.prof[i] = &seqdb_index.profile[i]->prof_sc;
  tree_align.align.path.compose (decomp, TRUE);
  return tree_align;
}

Progressive_aligner::~Progressive_aligner()
{ }

// Optimal_accuracy_progressive_aligner

Optimal_accuracy_progressive_aligner::Optimal_accuracy_progressive_aligner (Pair_PHMM& phmm, Dirichlet_prior& prior)
  : pscore (*prior.pscores),
    phmm (phmm),
    prior (prior),
    leaves (0),
    nodes (0)
{ }

Optimal_accuracy_progressive_aligner::~Optimal_accuracy_progressive_aligner()
{
  delete_post();
}

void Optimal_accuracy_progressive_aligner::delete_post()
{
  for (int x = 0; x < post.xsize(); ++x)
    for (int y = x + 1; y < post.ysize(); ++y)
      if (post(x,y))
	delete post(x,y);
}

void Optimal_accuracy_progressive_aligner::init_default_params()
{
  prior.initialise();
}

void Optimal_accuracy_progressive_aligner::load_params (const char* filename)
{
  Telegraph_PScores_adaptor tgio (prior);
  tgio.read (filename);
}

void Optimal_accuracy_progressive_aligner::train_params (const char* filename)
{
  // read in alignments
  Sequence_database stock_seq_db;
  Stockholm_database stock_db;
  ifstream stock_db_file (filename);
  if (!stock_db_file)
    THROWEXPR ("Pairwise alignment database '" << filename << "' not found");
  stock_db.read (stock_db_file, stock_seq_db);
  
  // digitize sequences
  stock_seq_db.seqs2dsqs (*phmm.alphabet);  // for counting symbol frequencies to get null model (see below)
  stock_seq_db.seqs2scores (*phmm.alphabet);  // for Pair_HMM DP
  
  // do EM
  Pair_HMM_trainer trainer (stock_db, phmm, prior, pscore);
  trainer.train();
}

void Optimal_accuracy_progressive_aligner::train_params_unaligned (const Sequence_database_index& seqdb_index)
{
  THROWEXPR ("train_params_unaligned() unimplemented");
}

void Optimal_accuracy_progressive_aligner::save_params (const char* filename)
{
  Telegraph_PScores_adaptor tgio (prior);
  tgio.write (filename);
}

void Optimal_accuracy_progressive_aligner::initialise (const Sequence_database_index& seqdb_index)
{
  // clear arrays
  post.resize (0, 0, (PostMat*) 0);
  expected_accuracy.resize (0, 0, 0.);

  // update sequence info
  leaves = seqdb_index.size();
  nodes = 2 * leaves - 1;
  leaf_name = seqdb_index.name;
  post.resize (leaves, leaves, (PostMat*) 0);
  expected_accuracy.resize (nodes, nodes, 0.);
  profile = seqdb_index.profile;

  leaf_seqlen = vector<int> (leaves);
  for (int i = 0; i < leaves; ++i)
    leaf_seqlen[i] = seqdb_index.profile[i]->size();

  // fill post
  fill_post();

  // do consistency transformation
  do_consistency_transformation();
}

void Optimal_accuracy_progressive_aligner::make_guide_tree()
{
  // reset tree & unaligned node vector
  tree = PHYLIP_tree();
  vector<int> unaligned (leaves);
  for (int i = 0; i < leaves; ++i)
    {
      unaligned[i] = i;
      tree.add_node();
    }

  // do greedy clustering
  while ((int) unaligned.size() > 1)
    {
      // print log message
      CTAG(2,OPTALIGN) << "Clusters: " << unaligned << "\n";

      // find best groups to merge
      double best_acc = -1.;
      int best_ui = -1;
      int best_uj = -1;
      for (int i = 0; i < (int) unaligned.size(); ++i)
	for (int j = i + 1; j < (int) unaligned.size(); ++j)
	  {
	    // get ui=unaligned[i] & uj=unaligned[j]
	    const int ui = unaligned[i];
	    const int uj = unaligned[j];
	    const double acc = expected_accuracy (ui, uj);
	    CTAG(1,OPTALIGN_MERGE) << "Expected accuracy of clusters " << ui << " and " << uj << " is " << acc << "\n";
	    if (acc > best_acc)
	      {
		best_ui = ui;
		best_uj = uj;
		best_acc = acc;
	      }
	  }

      // print log message
      CTAG(2,OPTALIGN) << "Merging clusters " << best_ui << " and " << best_uj << ", expected accuracy " << best_acc << "\n";

      // update tree
      const int a = tree.add_node();
      tree.add_branch (a, best_ui, 1.);
      tree.add_branch (a, best_uj, 1.);
      // remove best_ui & best_uj from unaligned
      vector<int> new_unaligned;
      new_unaligned.reserve (unaligned.size() - 1);
      for_const_contents (vector<int>, unaligned, ui)
	if (*ui != best_ui && *ui != best_uj)
	  new_unaligned.push_back (*ui);
      unaligned.swap (new_unaligned);
      // get expected_accuracy for new cluster
      for_const_contents (vector<int>, unaligned, uk)
	expected_accuracy (a, *uk) = expected_accuracy (*uk, a) = best_acc * (expected_accuracy (best_ui, *uk) + expected_accuracy (best_uj, *uk)) / 2;
      // add new cluster to unaligned
      unaligned.push_back (a);
    }

  // tidy up tree
  if (tree.nodes() != nodes)
    THROWEXPR ("Tree has " << tree.nodes() << " nodes, expected " << nodes);
  tree.root = nodes - 1;
  tree.rebuild_parents();
  tree.node_name = leaf_name;
  tree.node_name.insert (tree.node_name.end(), nodes - leaves, sstring());
}

Alignment_path Optimal_accuracy_progressive_aligner::align_axy (int a, int x, int y, double tx, double ty)
{
  // make vectors of leaf nodes descended from x & y, and Decomposition's for {x & descendants} and {y & descendants}
  Row_index_set xleaf, yleaf;
  Decomposition xdecomp, ydecomp;
  for_branches_post (tree, x, a, bi)
    {
      Row_pair rp (*bi);
      xdecomp[rp] = decomp[rp];
      if (tree.is_leaf ((*bi).second))
	xleaf.push_back ((*bi).second);
    }
  if (tree.is_leaf (x))
    xleaf.push_back (x);
  for_branches_post (tree, y, a, bi)
    {
      Row_pair rp (*bi);
      ydecomp[rp] = decomp[rp];
      if (tree.is_leaf ((*bi).second))
	yleaf.push_back ((*bi).second);
    }
  if (tree.is_leaf (y))
    yleaf.push_back (y);
  // make Alignment_path's for {x & descendants} and {y & descendants}
  Alignment_path xpath, ypath;
  if (xdecomp.size())
    xpath.compose (xdecomp, TRUE);
  else
    {
      xpath = Alignment_path (x+1, leaf_seqlen[x]);
      xpath[x] = Row (leaf_seqlen[x], TRUE);
    }
  if (ydecomp.size())
    ypath.compose (ydecomp, TRUE);
  else
    {
      ypath = Alignment_path (y+1, leaf_seqlen[y]);
      ypath[y] = Row (leaf_seqlen[y], TRUE);
    }
  // call align_subpaths()
  Pairwise_path xy_path = align_subpaths (xpath, xleaf, ypath, yleaf);
  // add an ancestral sequence with 1's in every column
  xy_path.insert_rows (0, Row (xy_path.columns(), 1));
  // return
  return xy_path;
}

void Optimal_accuracy_progressive_aligner::refine()
{
  // compose Alignment_path
  Alignment_path path;
  path.compose (decomp, TRUE);
  for (int refine_iter = 0; refine_iter < OPTACC_REFINEMENT_ITERATIONS; ++refine_iter)
    {
      // print log message
      CTAG(4,OPTALIGN) << "Starting refinement iteration #" << refine_iter+1 << "\n";

      // randomly partition leaf nodes into two vectors, x & y
      Row_index_set leaf (leaves), xleaf, yleaf;
      for (int i = 0; i < leaves; ++i)
	leaf[i] = i;
      for (int i = 0; i < leaves; ++i)
	{
	  const int j = i + Rnd::rnd_int (leaves - i);
	  if (i < leaves / 2)
	    xleaf.push_back (leaf[j]);
	  else
	    yleaf.push_back (leaf[j]);
	  leaf[j] = leaf[i];
	}

      // print leaf sets
      CTAG(3,OPTALIGN) << "Aligning subalignment of rows (" << xleaf << ") with subalignment of rows (" << yleaf << ")\n";
      
      // make Alignment_path's for x & y (leaf nodes only), erasing empty columns
      // x leaves...
      Alignment_path xpath (leaves, path.columns());
      for_const_contents (Row_index_set, xleaf, xl)
	xpath[*xl].swap (path[*xl]);
      xpath.erase_empty_columns();
      // y leaves...
      Alignment_path ypath (leaves, path.columns());
      for_const_contents (Row_index_set, yleaf, yl)
	ypath[*yl].swap (path[*yl]);
      ypath.erase_empty_columns();
      // call align_subpaths()
      Alignment_path new_path = align_subpaths (xpath, xleaf, ypath, yleaf);
      // combine x & y alignment paths
      vector<Row_index_set> rows;
      rows.push_back (xleaf);
      rows.push_back (yleaf);
      vector<Alignment_path> paths;
      paths.push_back (Subalignment_path (xpath, xleaf, FALSE));
      paths.push_back (Subalignment_path (ypath, yleaf, FALSE));
      new_path.explode (rows, paths);
      // insert empty rows for ancestors
      new_path.insert_rows (leaves, nodes - leaves);
      // fill ancestral rows
      for_rooted_nodes_post (tree, bi)
	{
	  // get children
	  const Node node = (*bi).second;
	  const Row_index_set children = tree.children (node, tree.parent[node]);  // assumes node == row, as always in this class...
	  if (children.size())
	    {
	      // add ancestral sequence with 1's in every column
	      for (int col = 0; col < new_path.columns(); ++col)
		for_const_contents (Row_index_set, children, c)
		  if (new_path (*c, col))
		    {
		      new_path[node][col] = TRUE;
		      break;
		    }
	    }
	}
      // swap path & new_path
      path.swap_path (new_path);
      // print log message
      if (CTAGGING(2,OPTALIGN))
	{
	  CL << "Alignment path after refinement:\n";
	  Alignment_path tmp_path (leaves);
	  for (int r = 0; r < leaves; ++r)  // swap tmp_path & leaf rows from path
	    tmp_path[r].swap (path[r]);
	  tmp_path.show (CL, leaf_name);
	  for (int r = 0; r < leaves; ++r)  // swap back
	    tmp_path[r].swap (path[r]);
	}
    }
  // update decomp
  for_rooted_nodes_post (tree, bi)
    {
      // get children
      const Node node = (*bi).second;
      const Row_index_set children = tree.children (node, tree.parent[node]);  // assumes node == row, as always in this class...
      // put pairwise node-to-child alignments in decomp
      for_const_contents (Row_index_set, children, c)
	decomp[Row_pair (node, *c)] = Pairwise_path (path, node, *c, TRUE);
    }
}

void Optimal_accuracy_progressive_aligner::fill_post()
{
  // evaluate HMM scores
  CTAG(2,OPTALIGN) << "Evaluating Pair HMM scores\n";
  Pair_HMM_scores hmm = phmm.eval_hmm_scores (pscore);

  // fill post
  for (int x = 0; x < leaves; ++x)
    for (int y = x + 1; y < leaves; ++y)
      {
	// make copies of vars
	const Named_profile& xnp = *profile[x];
	const Named_profile& ynp = *profile[y];
	const Score_profile& xprof = xnp.prof_sc;
	const Score_profile& yprof = ynp.prof_sc;
	const Biosequence& xseq = xnp.seq;
	const Biosequence& yseq = ynp.seq;
	const sstring& xname = xnp.name;
	const sstring& yname = ynp.name;

	// get posterior probabilities
	CTAG(2,OPTALIGN) << "Creating forward-backward matrix for '" << xname << "' vs '" << yname << "'\n";
	Pair_forward_backward_DP_matrix fb (hmm, xprof, yprof);
	CTAG(2,OPTALIGN) << "Creating posterior matrix for '" << xname << "' vs '" << yname << "'\n";
	post(x,y) = new Post_pair_HMM (fb);
	if (CTAGGING(1,OPTALIGN_DOTPLOT))
	  {
	    Pair_envelope pair_env;
	    pair_env.initialise_from_posterior_matrix (*post(x,y), 1./8.);  // use min prob of 1/8 for max ANSI color contrast
	    CTAG(1,OPTALIGN DOTPLOT) << "Posterior matrix for '" << xname << "' vs '" << yname << "'\n";  // refresh CTAG() flags, which may have been corrupted by Pair_envelope constructor
	    pair_env.render_dotplot (CL, xseq, yseq, 8);
	  }

	// find expected accuracy
	CTAG(2,OPTALIGN) << "Creating optimal accuracy matrix for '" << xname << "' vs '" << yname << "'\n";
	Optimal_accuracy_DP_matrix optacc (*post(x,y));
	if (CTAGGING(2,OPTALIGN))
	  {
	    Pairwise_path xy_path = optacc.traceback();
	    Alignment xy_align (xy_path, xnp, ynp);
	    xy_align.write_MUL (CL, *phmm.alphabet);
	  }
	// store expected accuracy
	const double acc = optacc.final_score();
	expected_accuracy(x,y) = expected_accuracy(y,x) = acc;
	CTAG(4,OPTALIGN) << "Expected accuracy for '" << xname << "' vs '" << yname << "' is " << acc << "\n";
      }
}


void Optimal_accuracy_progressive_aligner::do_consistency_transformation()
{
  // do consistency transformation
  for (int n = 0; n < OPTACC_CONSISTENCY_ITERATIONS; ++n)
    {
      CTAG(4,OPTALIGN) << "Starting consistency transformation #" << n+1 << "\n";
      array2d<PostMat*> new_post (leaves, leaves, (Post_pair_HMM*) 0);
      for (int x = 0; x < leaves; ++x)
	for (int y = x + 1; y < leaves; ++y)
	  {
	    // make copies of vars
	    const Named_profile& xnp = *profile[x];
	    const Named_profile& ynp = *profile[y];
	    const int xsize = leaf_seqlen[x];
	    const int ysize = leaf_seqlen[y];
	    const sstring& xname = xnp.name;
	    const sstring& yname = ynp.name;

	    // print log message
	    CTAG(3,OPTALIGN) << "Starting consistency transformation #" << n+1 << " for '" << xname << "' vs '" << yname << "'\n";

	    // create transformed reward matrix
	    new_post (x, y) = new array2d<Prob> (xsize, ysize, 0.);
	    PostMat& new_xy (*new_post (x, y));

	    // handle z < x
	    for (int z = 0; z < x; ++z)
	      {
		const int zsize = leaf_seqlen[z];
		PostMat& old_zx (*post (z, x));
		PostMat& old_zy (*post (z, y));
		for (int i = 0; i < xsize; ++i)
		  for (int j = 0; j < ysize; ++j)
		    for (int k = 0; k < zsize; ++k)
		      new_xy(i,j) += old_zx(k,i) * old_zy(k,j);
	      }

	    // handle x < z < y
	    for (int z = x + 1; z < y; ++z)
	      {
		const int zsize = leaf_seqlen[z];
		PostMat& old_xz (*post (x, z));
		PostMat& old_zy (*post (z, y));
		for (int i = 0; i < xsize; ++i)
		  for (int j = 0; j < ysize; ++j)
		    for (int k = 0; k < zsize; ++k)
		      new_xy(i,j) += old_xz(i,k) * old_zy(k,j);
	      }

	    // handle z > y
	    for (int z = y + 1; z < leaves; ++z)
	      {
		const int zsize = leaf_seqlen[z];
		PostMat& old_xz (*post (x, z));
		PostMat& old_yz (*post (y, z));
		for (int i = 0; i < xsize; ++i)
		  for (int j = 0; j < ysize; ++j)
		    for (int k = 0; k < zsize; ++k)
		      new_xy(i,j) += old_xz(i,k) * old_yz(j,k);
	      }

	    // divide transformed posterior probabilities by (number of leaves minus one)
	    for (int i = 0; i < xsize; ++i)
	      for (int j = 0; j < ysize; ++j)
		new_xy(i,j) /= (double) (leaves - 1);
	  }
      // update post
      delete_post();
      post.swap (new_post);
    }
}

Pairwise_path Optimal_accuracy_progressive_aligner::align_subpaths (const Alignment_path& xpath, const Row_index_set& xrows,
								    const Alignment_path& ypath, const Row_index_set& yrows)
{
  // print log message
  if (CTAGGING(2,OPTALIGN))
    {
      // create subpaths
      const Subalignment_path xsubpath (xpath, xrows, TRUE);
      const Subalignment_path ysubpath (ypath, yrows, TRUE);
      // get row names
      vector<sstring> xname, yname;
      xname.reserve (xrows.size());
      for_const_contents (Row_index_set, xrows, xrow)
	xname.push_back (leaf_name[*xrow]);
      yname.reserve (yrows.size());
      for_const_contents (Row_index_set, yrows, yrow)
	yname.push_back (leaf_name[*yrow]);
      // print
      CTAG(2,OPTALIGN) << "Aligning subalignments X (rows " << xrows << ") and Y (rows " << yrows << ")\n";  // reset log tags
      CL << "Subpath X:\n";
      xsubpath.show (CL, xname);
      CL << "Subpath Y:\n";
      ysubpath.show (CL, yname);
    }

  // create the expected-sum-of-pairs-score reward matrix
  array2d<Prob> reward (xpath.columns(), ypath.columns(), 0.);
  vector<int> xcoords = xpath.create_seq_coords();
  for (int x = 0; x < xpath.columns(); ++x)
    {
      vector<int> ycoords = ypath.create_seq_coords();
      for (int y = 0; y < ypath.columns(); ++y)
	{
	  Prob sps = 0.;
	  for_const_contents (Row_index_set, xrows, xrow)
	    for_const_contents (Row_index_set, yrows, yrow)
	    if (*xrow < *yrow)
	      sps += (*post(*xrow,*yrow)) (xcoords[*xrow], ycoords[*yrow]);
	    else
	      sps += (*post(*yrow,*xrow)) (ycoords[*yrow], xcoords[*xrow]);
	  reward (x, y) = sps;
	  ypath.inc_seq_coords (ycoords, y);
	}
      xpath.inc_seq_coords (xcoords, x);
    }
  // build optimal accuracy DP matrix
  Optimal_accuracy_DP_matrix optacc (reward);
  // return traceback
  return optacc.traceback();
}
