#ifndef COVMODEL_INCLUDED
#define COVMODEL_INCLUDED

#include "scfg/pairpcfg.h"

struct CM_node : Pair_CFG_state_type_enum, Grammar_state_enum
{
  // data
  int base_state;
  vector<CM_node*> child;
  // destructor
  virtual ~CM_node() { }
  // state count & types
  virtual int states() const = 0;
  virtual State_type state_type (int state) const = 0;
  virtual State_type pairwise_state_type (int state) const = 0;
  virtual int consensus_emit_len() const { return 0; }
  // parameter builder methods
  virtual void build_params (PScope& pscope) = 0;
  virtual void build_prior (Dirichlet_prior& prior, const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul) = 0;
  // CFG builder methods
  // connect() calls dest_node.connect_incoming_transitions_from_*()
  virtual void connect_incoming_transitions_from_insert (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const = 0;
  virtual void connect_incoming_transitions_from_delete (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const = 0;
  virtual void connect_incoming_transitions_from_match (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const = 0;
  virtual void connect (Pair_PCFG& cfg) const = 0;
  // connect_single & connect_pairwise methods call cfg.init_emit, followed by connect()
  void connect_single (Pair_PCFG& cfg) const;
  void connect_pairwise (Pair_PCFG& cfg) const;
  // Node name
  virtual const char* cm_node_type() = 0;
};

struct CM_end : CM_node
{
  // CM_node methods
  int states() const { return 0; }
  State_type state_type (int state) const { return Undefined; }
  State_type pairwise_state_type (int state) const { return Undefined; }
  void build_params (PScope& pscope) { }
  void build_prior (Dirichlet_prior& prior, const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul) { }
  void connect_incoming_transitions_from_insert (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    connect_incoming_transitions_from_match (cfg, src_state, src_prob);
  }
  void connect_incoming_transitions_from_delete (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    connect_incoming_transitions_from_match (cfg, src_state, src_prob);
  }
  void connect_incoming_transitions_from_match (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    cfg.transition (src_state, End) = src_prob;
  }
  void connect (Pair_PCFG& cfg) const { }
  // node type descriptor
  const char* cm_node_type() { return Grammar_end_state_name; }
};

struct CM_bifurc : CM_node
{
  // CM_node methods
  int states() const { return 3; }

  State_type state_type (int state) const
  {
    switch (state - base_state)
      {
      case 0:  // b
	return Bifurc;
      case 1:  // l
      case 2:  // r
	return Null;
      default: break;
      }
    return Undefined;
  };

  State_type pairwise_state_type (int state) const
  {
    return state_type (state);
  };

  void build_params (PScope& pscope) { }
  void build_prior (Dirichlet_prior& prior, const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul) { }
  void connect_incoming_transitions_from_insert (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const { }
  void connect_incoming_transitions_from_delete (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const { }
  void connect_incoming_transitions_from_match (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const { }
  void connect (Pair_PCFG& cfg) const
  {
    cfg.init_bifurc (b(), l(), r());
    const PFunc one (1.);
    child[0]->connect_incoming_transitions_from_match (cfg, l(), one);
    child[1]->connect_incoming_transitions_from_match (cfg, r(), one);
  }
  // CM_bifurc state index accessor
  int b() const { return base_state; }
  int l() const { return base_state + 1; }
  int r() const { return base_state + 2; }
  // node type descriptor
  const char* cm_node_type() { return "Bifurcation"; }
};

struct CM_start : CM_node
{
  // CM_node methods
  int states() const { return 0; }

  State_type state_type (int state) const { return Undefined; }
  State_type pairwise_state_type (int state) const { return Undefined; }

  void build_params (PScope& pscope) { }
  void build_prior (Dirichlet_prior& prior, const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul) { }
  void connect_incoming_transitions_from_insert (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const { }
  void connect_incoming_transitions_from_delete (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const { }
  void connect_incoming_transitions_from_match (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const { }
  void connect (Pair_PCFG& cfg) const
  {
    child[0]->connect_incoming_transitions_from_match (cfg, start_state(), PFunc(1.));
  }
  // CM_start state index accessor
  int start_state() const { return Start; }
  // node type descriptor
  const char* cm_node_type() { return "Start"; }
};

struct CM_emit_single : CM_node
{
  // CM_emit_single data
  PGroup m2md, d2md, m2ix, i2ix;
  Alphabet_group mat, ins;

  // CM_node methods
  int states() const { return 3; }

  State_type state_type (int state) const
  {
    switch (state - base_state)
      {
      case 0:  // m
      case 1:  // i
	return x_emit_state_type();
      case 2:  // d
	return Null;
      default: break;
      }
    return Undefined;
  };

  State_type pairwise_state_type (int state) const
  {
    switch (state - base_state)
      {
      case 0:  // m
	return xy_emit_state_type();
      case 1:  // i
	return y_emit_state_type();
      case 2:  // d
	return x_emit_state_type();
      default: break;
      }
    return Undefined;
  };

  int consensus_emit_len() const { return 1; }

  void build_params (PScope& pscope)
  {
    m2md = pscope.new_group (2, "m2md", "md");
    d2md = pscope.new_group (2, "d2md", "md");
    m2ix = pscope.new_group (2, "m2ix", "ix");
    i2ix = pscope.new_group (2, "i2ix", "ix");
    mat = pscope.new_alphabet_group (CFG_alphabet, "mat");
    ins = pscope.new_alphabet_group (CFG_alphabet, "ins");
  }
  void build_prior (Dirichlet_prior& prior, const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul)
  {
    prior.assign_Laplace (m2md);
    prior.assign_Laplace (d2md);
    prior.assign_Laplace (m2ix);
    prior.assign_Laplace (i2ix);

    Dirichlet_mixture single_dm (mat.group_size, 1);
    for (int c = 0; c < mat.group_size; ++c)
      single_dm.alpha[0][0][c] = mul * FScore2Prob (ps[single_emit[c]]);
    prior.assign (mat, single_dm);
    prior.assign (ins, single_dm);
  }
  void connect_incoming_transitions_from_insert (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    connect_incoming_transitions_from_match (cfg, src_state, src_prob);
  }
  void connect_incoming_transitions_from_match (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    cfg.transition (src_state, m()) = src_prob * m2md[0];
    cfg.transition (src_state, d()) = src_prob * m2md[1];
  }
  void connect_incoming_transitions_from_delete (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    cfg.transition (src_state, m()) = src_prob * d2md[0];
    cfg.transition (src_state, d()) = src_prob * d2md[1];
  }
  void connect (Pair_PCFG& cfg) const
  {
    for (int c = 0; c < CFG_alphabet_size; ++c)
      {
	cfg.emit[m()][c] = mat[c];
	cfg.emit[i()][c] = ins[c];
      }

    cfg.transition (m(), i()) = m2ix[0];
    cfg.transition (i(), i()) = i2ix[0];

    CM_node& dest_node (*child[0]);
    dest_node.connect_incoming_transitions_from_match (cfg, m(), m2ix[1]);
    dest_node.connect_incoming_transitions_from_insert (cfg, i(), i2ix[1]);
    dest_node.connect_incoming_transitions_from_delete (cfg, d(), PFunc(1.));
  }
  // CM_emit_single state index accessors
  int m() const { return base_state; }
  int i() const { return base_state + 1; }
  int d() const { return base_state + 2; }

  // virtual CM_emit_single emit state type methods
  virtual State_type x_emit_state_type() const = 0;
  virtual State_type y_emit_state_type() const = 0;
  virtual State_type xy_emit_state_type() const = 0;
};

struct CM_emit_l : CM_emit_single
{
  State_type x_emit_state_type() const { return EmitXL; }
  State_type y_emit_state_type() const { return EmitYL; }
  State_type xy_emit_state_type() const { return EmitXLYL; }
  // node type descriptor
  const char* cm_node_type() { return "EmitL"; }
};

struct CM_emit_r : CM_emit_single
{
  State_type x_emit_state_type() const { return EmitXR; }
  State_type y_emit_state_type() const { return EmitYR; }
  State_type xy_emit_state_type() const { return EmitXRYR; }
  // node type descriptor
  const char* cm_node_type() { return "EmitR"; }
};

struct CM_emit_pair : CM_node
{
  // CM_emit_pair data
  PGroup m2md, d2md, m2ijkx, i2ijkx, j2jkx, k2kx;
  Alphabet_group mat, insl, insr, inslr;

  // CM_node methods
  int states() const { return 5; }

  State_type state_type (int state) const
  {
    switch (state - base_state)
      {
      case 0:  // m
      case 1:  // i
	return EmitXLR;
      case 2:  // j
	return EmitXL;
      case 3:  // k
	return EmitXR;
      case 4:  // d
	return Null;
      default: break;
      }
    return Undefined;
  };

  State_type pairwise_state_type (int state) const
  {
    switch (state - base_state)
      {
      case 0:  // m
	return EmitXLRYLR;
      case 1:  // i
	return EmitYLR;
      case 2:  // j
	return EmitYL;
      case 3:  // k
	return EmitYR;
      case 4:  // d
	return EmitXLR;
      default: break;
      }
    return Undefined;
  };

  int consensus_emit_len() const { return 2; }

  void build_params (PScope& pscope)
  {
    m2md = pscope.new_group (2, "m2md", "md");
    d2md = pscope.new_group (2, "d2md", "md");
    m2ijkx = pscope.new_group (4, "m2ijkx", "ijkx");
    i2ijkx = pscope.new_group (4, "i2ijkx", "ijkx");
    j2jkx = pscope.new_group (3, "j2jkx", "jkx");
    k2kx = pscope.new_group (2, "k2kx", "kx");
    mat = pscope.new_alphabet_group (CFG_alphabet, 2, "matlr");
    insl = pscope.new_alphabet_group (CFG_alphabet, "insl");
    insr = pscope.new_alphabet_group (CFG_alphabet, "insr");
    inslr = pscope.new_alphabet_group (CFG_alphabet, 2, "inslr");
  }
  void build_prior (Dirichlet_prior& prior, const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul)
  {
    prior.assign_Laplace (m2md);
    prior.assign_Laplace (d2md);
    prior.assign_Laplace (m2ijkx);
    prior.assign_Laplace (i2ijkx);
    prior.assign_Laplace (j2jkx);
    prior.assign_Laplace (k2kx);

    Dirichlet_mixture pair_dm (mat.group_size, 1);
    Dirichlet_mixture single_dm (insl.group_size, 1);
    for (int c = 0; c < insl.group_size; ++c)
      single_dm.alpha[0][0][c] = mul * FScore2Prob (ps[single_emit[c]]);
    for (int c = 0; c < mat.group_size; ++c)
      pair_dm.alpha[0][0][c] = mul * FScore2Prob (ps[pair_emit[c]]);
    prior.assign (mat, pair_dm);
    prior.assign (insl, single_dm);
    prior.assign (insr, single_dm);
    prior.assign (inslr, pair_dm);
  }
  void connect_incoming_transitions_from_insert (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    connect_incoming_transitions_from_match (cfg, src_state, src_prob);
  }
  void connect_incoming_transitions_from_match (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    cfg.transition (src_state, m()) = src_prob * m2md[0];
    cfg.transition (src_state, d()) = src_prob * m2md[1];
  }
  void connect_incoming_transitions_from_delete (Pair_PCFG& cfg, int src_state, const PFunc& src_prob) const
  {
    cfg.transition (src_state, m()) = src_prob * d2md[0];
    cfg.transition (src_state, d()) = src_prob * d2md[1];
  }
  void connect (Pair_PCFG& cfg) const
  {
    vector<int> c (2);
    int& cl (c[0]);
    int& cr (c[1]);
    const int lmul = cfg.emit_xl_mul (state_type(m()));
    const int rmul = cfg.emit_xr_mul (state_type(m()));
    for (cl = 0; cl < CFG_alphabet_size; ++cl)
      {
	cfg.emit[j()][cl] = insl[cl];
	cfg.emit[k()][cl] = insr[cl];
	for (cr = 0; cr < CFG_alphabet_size; ++cr)
	  {
	    const int idx = mat.intvec2index (c);
	    cfg.emit[m()][cl*lmul + cr*rmul] = mat[idx];
	    cfg.emit[i()][cl*lmul + cr*rmul] = inslr[idx];
	  }
      }

    cfg.transition (m(), i()) = m2ijkx[0];
    cfg.transition (m(), j()) = m2ijkx[1];
    cfg.transition (m(), k()) = m2ijkx[2];

    cfg.transition (i(), i()) = i2ijkx[0];
    cfg.transition (i(), j()) = i2ijkx[1];
    cfg.transition (i(), k()) = i2ijkx[2];

    cfg.transition (j(), j()) = j2jkx[0];
    cfg.transition (j(), k()) = j2jkx[1];

    cfg.transition (k(), k()) = k2kx[0];

    CM_node& dest_node (*child[0]);
    dest_node.connect_incoming_transitions_from_match (cfg, m(), m2ijkx[3]);
    dest_node.connect_incoming_transitions_from_insert (cfg, i(), i2ijkx[3]);
    dest_node.connect_incoming_transitions_from_insert (cfg, j(), j2jkx[2]);
    dest_node.connect_incoming_transitions_from_insert (cfg, k(), k2kx[1]);
    dest_node.connect_incoming_transitions_from_delete (cfg, d(), PFunc(1.));
  }
  // CM_emit_pair state index accessors
  int m() const { return base_state; }
  int i() const { return base_state + 1; }
  int j() const { return base_state + 2; }
  int k() const { return base_state + 3; }
  int d() const { return base_state + 4; }
  // node type descriptor
  const char* cm_node_type() { return "EmitLR"; }
};

struct Covariance_model : Pair_CFG_state_type_enum
{
  // data
  vector<CM_node*> node;
  vector<int> state2node;
  PScores pscores;
  Dirichlet_prior prior;

  // models, consensus seq
  Named_profile consensus_seq;
  Pair_PCFG single_model;
  Pair_PCFG pair_model;  // this exists only for the state typing; the PFunc's herein are meaningless

  // constructor
  Covariance_model();
  // destructor
  ~Covariance_model();

  // builder methods
  void build_from_parse_tree (const vector<State_type>& parse_tree);
  void init_prior (const PGroup& pair_emit, const PGroup& single_emit, const PScores& ps, double mul);
  void add_node (CM_node* new_node, stack<CM_node*>& parent, bool bifurc);

  // helpers
  int states() const;

  // the following method converts a local traceback for an emitted sequence
  // into a pairwise alignment to a "consensus" sequence for the model.
  // Row 0 of the alignment represents the consensus; row 1 represents the emitted sequence.
  Pairwise_path convert_local_path_to_consensus_alignment (const Pair_CFG_local_path& local_path, const Named_profile& np);

  // show method
  void show (ostream& out) const;
};

#endif /* COVMODEL_INCLUDED */
