#include "hmm/singledp.h"
#include "util/dexception.h"

void _throw_unimpl()
{
  THROW Standard_exception ("Method unimplemented");
}

int Single_Viterbi_interface::get_Viterbi_score() const { _throw_unimpl(); return -InfinityScore; }
vector<int> Single_Viterbi_interface::optimal_state_path() const { _throw_unimpl(); return vector<int>(); }

int Single_forward_interface::get_forward_score() const { _throw_unimpl(); return -InfinityScore; }
vector<int> Single_forward_interface::sample_state_path() const { _throw_unimpl(); return vector<int>(); }

Score Single_forward_backward_interface::get_forward_score() const { _throw_unimpl(); return -InfinityScore; }
const Single_HMM_counts& Single_forward_backward_interface::get_expected_counts() const { _throw_unimpl(); return *((const Single_HMM_counts*) 0); }
const vector<Metaprob>& Single_forward_backward_interface::get_expected_metacounts() const { _throw_unimpl(); return *((const vector<Metaprob>*) 0); }

Single_Viterbi_interface::~Single_Viterbi_interface() { }

Single_forward_interface::~Single_forward_interface() { }

Single_forward_backward_interface::~Single_forward_backward_interface() { }

// method to display a set of metascore vectors
void Single_forward_backward_interface::show_metacounts (ostream& o) const
{
  int old_prec = o.precision(4);
  save_flags (o);
  right_align (o);

  const vector<Metaprob>& metacounts = get_expected_metacounts();
  if (metacounts.size())
    {
      o << "Metacounts:\n";
      o << "Sequence pos:";
      for (int j = 0; j < (int) metacounts[0].size(); ++j)
	{
	  o << " ";
	  o.width(10);
	  o << j;
	}
      o << "\n";
      for (int idx = 0; idx < (int) metacounts.size(); ++idx)
	{
	  o << "Metaindex ";
	  o.width(2);
	  o << idx << ":";
	  for (int j = 0; j < (int) metacounts[idx].size(); ++j)
	    {
	      o << " ";
	      o.width(10);
	      o << metacounts[idx][j];
	    }
	  o << "\n";
	}
    }

  restore_flags (o);
  o.precision (old_prec);
}

