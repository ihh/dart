// pk 3/05 accept long input lines in stockholm file

#include "util/logfile.h"
#include "util/rnd.h"
#include "irrev/irrev_em_matrix.h"
#include "util/vector_output.h"

// main program
int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  INIT_OPTS_LIST (opts, argc, argv, 2, "[options]  <rate matrix filename>  <tree filename>",
		  "simulate a multiple alignment from a substitution model, keeping track of EM statistics\n");

  opts.print_title ("Simulation options");

  int cols;
  opts.add ("c -columns", cols = 100, "number of alignment columns", FALSE);

  opts.newline();
  opts.print_title ("File formats");
  opts.print ("Rate matrix is in xrate format; tree is in New Hampshire format.\n");
  opts.print ("Ancestral sequences at internal nodes of the tree will not be printed,\n");
  opts.print ("*unless* they are labeled with names in the New Hampshire tree file.\n");

  // parse the command line
  //
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  try
    {
      // get args
      const char* irrev_filename = opts.args[0].c_str();
      const char* tree_filename = opts.args[1].c_str();

      // read rate matrix
      Irrev_EM_matrix irrev (1, 1);
      ifstream irrev_file (irrev_filename);
      if (!irrev_file)
	THROWEXPR ("Rate matrix file '" << irrev_filename << "' not found");
      irrev.read (irrev_file);

      // read tree
      PHYLIP_tree tree;
      ifstream tree_file (tree_filename);
      if (!tree_file)
	THROWEXPR ("Tree file '" << tree_filename << "' not found");
      tree.read (tree_file);

      // get alphabet
      const Alphabet& alph = irrev.alphabet();

      // figure out which nodes have names
      int named_nodes = 0;
      vector<int> node2row (tree.nodes(), -1);
      for (int node = 0; node < tree.nodes(); ++node)
	if (tree.node_name[node].size())
	  node2row[node] = named_nodes++;

      // initialise Stockade
      Stockade stock (named_nodes, cols);
      for (int node = 0; node < tree.nodes(); ++node)
	if (node2row[node] >= 0)
	  stock.align.row_name[node2row[node]] = tree.node_name[node];
      for (int row = 0; row < named_nodes; ++row)
	for (int col = 0; col < cols; ++col)
	  stock.align.path.row(row)[col] = 1;

      // initialise Update_statistics
      const int total_states = irrev.m();
      Update_statistics stats (total_states);

      // get rate matrix
      const Matrix& R (irrev.R);

      // loop over alignment columns
      for (int col = 0; col < cols; ++col)
	{
	  // print log message
	  CTAG(6,XSIM) << "Simulating column " << col << "\n";

	  // keep track of state at each node
	  vector<int> node_state (tree.nodes());

	  // sample state of root node
	  const int root_state = Rnd::choose (irrev.pi);
	  node_state[tree.root] = root_state;
	  ++stats.s[root_state];
	  if (node2row[tree.root] >= 0)
	    stock.np[node2row[tree.root]].seq.push_back (alph.int2char (root_state));

	  // print log message
	  CTAG(5,XSIM) << "Root state (node " << tree.root << ") is " << root_state << "\n";

	  // simulate evolution along each branch
	  for_rooted_branches_pre (tree, b)
	    {
	      const int parent = (*b).first;
	      const int child = (*b).second;
	      double time_remaining = (*b).length;

	      // get state of parent node
	      int state = node_state[parent];

	      // print log message
	      CTAG(5,XSIM) << "Simulating branch from node " << parent << " to " << child
			   << ", time " << time_remaining << ", starting in state " << state << "\n";

	      // sample a trajectory for the Markov chain along the branch
	      while (true)
		{
		  // get total rate
		  const double total_rate = -R(state+1,state+1);
		  CTAG(3,XSIM) << "Total mutation rate from state " << state << " is " << total_rate << "\n";

		  // sample time to next mutation
		  const double rnd = Rnd::prob();
		  const double time_to_next = -Math_fn::math_log (rnd) / total_rate;

		  // did a mutation occur before the end of the branch?
		  if (time_remaining > time_to_next)
		    {
		      // get row of rate matrix
		      vector<double> row (total_states, 0.);
		      for (int i = 0; i < total_states; ++i)
			if (i != state)
			  row[i] = R(state+1,i+1);

		      // print log message
		      CTAG(3,XSIM) << "Row " << state << " of rate matrix is (" << row << ")\n";

		      // sample next state
		      const int next_state = Rnd::choose (row);

		      // update stats
		      ++stats.u (state, next_state);
		      stats.w[state] += time_to_next;

		      // print log message
		      CTAG(4,XSIM) << "Mutation from " << state << " to " << next_state
				   << " after wait " << time_to_next << " (time remaining " << time_remaining - time_to_next << ")\n";

		      // update state, time
		      state = next_state;
		      time_remaining -= time_to_next;
		    }
		  else // no mutation
		    {
		      // print log message
		      CTAG(4,XSIM) << "No mutation in remaining time\n";
		  
		      // update stats
		      stats.w[state] += time_remaining;

		      // exit loop
		      break;
		    }
		}

	      // set state of child node
	      node_state[child] = state;
	      if (node2row[child] >= 0)
		stock.np[node2row[child]].seq.push_back (alph.int2char (state));

	      // print log message
	      CTAG(5,XSIM) << "Child node " << child << " finished in state " << state << "\n";
	    }
	}

      // add tree to Stockholm alignment
      sstring tree_string;
      tree.write (tree_string, 0);
      tree_string.chomp();
      stock.align.add_gf_annot (sstring (Stockholm_New_Hampshire_tag), tree_string);

      // add Update_statistics to Stockholm alignment
      sstring stats_string;
      stats_string << stats;  // dump Update_statistics to string
      vector<sstring> stats_lines = stats_string.split ("\n");
      stats_lines.erase (stats_lines.begin());   // HACK: shave off first line of Update_statistics dump string, which contains irrelevant log-likelihood
      for_const_contents (vector<sstring>, stats_lines, stats_line)
	stock.align.add_gf_annot (sstring (Stockholm_comment_tag), *stats_line);

      // output alignment
      stock.align.write_Stockholm (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << "ERROR: " << e.what();
      exit(1);
    }

  return 0;
}
