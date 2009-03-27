#include "scfg/covmodel.h"

#define CM_CONSENSUS_SEQ_NAME "Consensus"

void CM_node::connect_single (Pair_PCFG& cfg) const
{
  const PFunc one (1.);
  for (int s = base_state; s < base_state + states(); ++s)
    cfg.init_emit (s, state_type(s), one);
  connect (cfg);
}

void CM_node::connect_pairwise (Pair_PCFG& cfg) const
{
  const PFunc one (1.);
  for (int s = base_state; s < base_state + states(); ++s)
    cfg.init_emit (s, pairwise_state_type(s), one);
  connect (cfg);
}

Covariance_model::Covariance_model() : pscores(), prior (pscores), single_model (0), pair_model (0)
{ }

Covariance_model::~Covariance_model()
{
  for_contents (vector<CM_node*>, node, n)
    delete *n;
}

void Covariance_model::build_from_parse_tree (const vector<State_type>& parse_tree)
{
  stack<CM_node*> parent;
  add_node (new CM_start(), parent, false);
  for_const_contents (vector<State_type>, parse_tree, t)
    {
      switch (*t)
	{
	case EmitXLYL:
	case EmitXL:
	case EmitYL:
	  add_node (new CM_emit_l(), parent, false);
	  break;
	case EmitXRYR:
	case EmitXR:
	case EmitYR:
	  add_node (new CM_emit_r(), parent, false);
	  break;
	case EmitXLRYLR:
	case EmitXLR:
	case EmitYLR:
	case EmitXLYLR:
	case EmitXRYLR:
	case EmitXLRYL:
	case EmitXLRYR:
	case EmitXLYR:
	case EmitXRYL:
	  add_node (new CM_emit_pair(), parent, false);
	  break;
	case Bifurc:
	  add_node (new CM_bifurc(), parent, true);
	  break;
	case Null:
	  add_node (new CM_end(), parent, false);
	  parent.pop();
	  break;
	default:
	  break;
	}
    }

  if (!parent.empty())
    THROWEXPR ("parent stack not empty");

  // build models
  single_model = Pair_PCFG (states());
  for_const_contents (vector<CM_node*>, node, n)
    (*n)->connect_single (single_model);

  pair_model = Pair_PCFG (states());
  for_const_contents (vector<CM_node*>, node, n)
    (*n)->connect_pairwise (pair_model);

  // get consensus seq len
  int consensus_seq_len = 0;
  for_const_contents (vector<CM_node*>, node, n)
    consensus_seq_len += (*n)->consensus_emit_len();

  // build consensus sequence
  // some parts of Named_profile left uninitialized here...
  consensus_seq.clear();
  consensus_seq.name = CM_CONSENSUS_SEQ_NAME;
  consensus_seq.seq = sstring (consensus_seq_len, (char) '*');
  consensus_seq.dsq = Digitized_biosequence (consensus_seq_len);
  consensus_seq.prof_sc = Score_profile (consensus_seq_len, CFG_alphabet.flat_score_map (0));
}

void Covariance_model::init_prior (const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul)
{
  for_const_contents (vector<CM_node*>, node, n)
    (*n)->build_prior (prior, pair_emit, single_emit, ps, mul);
}

void Covariance_model::add_node (CM_node* new_node, stack<CM_node*>& parent, bool bifurc)
{
  node.push_back (new_node);
  new_node->base_state = states();
  const int new_states = new_node->states();
  for (int i = 0; i < new_states; ++i)
    state2node.push_back (node.size() - 1);
  new_node->build_params (pscores);
  if (!parent.empty())
    {
      parent.top()->child.push_back (new_node);
      parent.pop();
    }
  parent.push (new_node);
  if (bifurc)
    parent.push (new_node);
}

int Covariance_model::states() const
{
  return state2node.size();
}

Pairwise_path Covariance_model::convert_local_path_to_consensus_alignment (const Pair_CFG_local_path& local_path, const Named_profile& np)
{
  // this is a very messily implemented function.... then again, it's doing a fairly messy job.
  // Input is a local alignment of a sequence to a covariance model, together with the sequence itself (np).
  // Output is a pairwise alignment between the "consensus" sequence for the covariance model and np.

  // convert single local path to pairwise local path
  Pair_CFG_local_path pair_local_path;
  pair_local_path.xstart = 0;
  pair_local_path.xlen = consensus_seq.size();
  pair_local_path.ystart = local_path.xstart;
  pair_local_path.ylen = local_path.xlen;
  pair_local_path.path = local_path.path;
  
  // create a pair SCFG
  const Pair_CFG_scores pair_cfg = pair_model.eval_cfg_scores (pscores);
  if (CTAGGING(1,COVMODEL)) {
    CL << "Covariance model Pair_CFG_scores:\n";
    pair_cfg.show (CL);
    CL << "Pair_CFG_local_path:\n";
    pair_local_path.show (CL);
  }

  // use the pair SCFG to turn the pairwise local path into a consensus-to-sequence alignment
  const Pair_CFG_parse_tree pair_parse = pair_cfg.parse (pair_local_path);
  const Pair_CFG_alignment pair_align = pair_parse.alignment (pair_cfg.state_type, consensus_seq, np);

  // return alignment
  return pair_align.pairwise_path();
}

void Covariance_model::show (ostream& out) const
{
  out << "Covariance model nodes:\n";
  for (int i = 0; i < (int) node.size(); ++i)
    out << " Node " << i << " (" << node[i]->cm_node_type() << ") Base state " << node[i]->base_state << "\n";
  out << "Consensus sequence length = " << consensus_seq.size() << "\n";
  out << "Covariance model single SCFG\n";
  single_model.show (out);
  out << "Covariance model pair SCFG\n";
  pair_model.show (out);
}
