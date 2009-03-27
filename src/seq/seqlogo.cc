#include "seq/seqlogo.h"

Sequence_logo::Sequence_logo (const Score_profile& prof_sc, const Alphabet& alphabet)
  : alphabet(alphabet),
    text_height(16)
{
  Weight_profile prof_w = prof_sc.score2weight();
  prof_w.normalise();
  for (int col = 0; col < (int) prof_w.size(); ++col)
    {
      sort_prof_w.push_back (vector<Symbol_weight> (prof_w[col].begin(), prof_w[col].end()));
      vector<Symbol_weight>& col_data = sort_prof_w[col];
      // the following is really ugly
      vector<double> sym_wt (col_data.size());
      vector<int> sym_order (col_data.size());
      for (int i = 0; i < (int) col_data.size(); ++i)
	{
	  sym_wt[i] = col_data[i].second;
	  sym_order[i] = i;
	}
      Schwartzian<double> by_sym_wt (sym_wt);
      sort (sym_order.begin(), sym_order.end(), by_sym_wt);
      vector<Symbol_weight> tmp_col_data (col_data);
      for (int i = 0; i < (int) sym_order.size(); ++i)
	col_data[i] = tmp_col_data[sym_order[i]];
    }
}

void Sequence_logo::show (ostream& o) const
{
  for (int row = 0; row < text_height; ++row)
    {
      const double p = ((double) row) / ((double) text_height);
      if (row % (text_height / 4) == 0)
	{
	  if (p == 0) o << "100 -- ";
	  else if (p == .25) o << " 75 -  ";
	  else if (p == .5) o << " 50 -- ";
	  else if (p == .75) o << " 25 -  ";
	}
      else o << "       ";
      for (int col = 0; col < (int) sort_prof_w.size(); ++col)
	{
	  double q = p;
	  int i = 0;
	  while (i < (int) sort_prof_w[col].size() - 1)
	    if ((q -= sort_prof_w[col][i].second) > 0) break;
	    else ++i;
	  o << alphabet.int2char (sort_prof_w[col][i].first);
	}
      o << '\n';
    }
}
