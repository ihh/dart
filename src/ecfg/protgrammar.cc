#include "ecfg/protgrammar.h"
#include "tree/pam.h"

Protein_grammar::Protein_grammar (int zones, int hidden)
  : ECFG_scores (Protein_alphabet, zones)
{
  PAM_factory pam_factory;
  const array2d<double> pam_matrix = pam_factory.rate_matrix();
  const vector<double> pam_prior = pam_factory.prior();

  for (int s = 0; s < zones; ++s)
    {
      // create a dummy EM_matrix just to get the alphabet (eww)
      EM_matrix dummy_em_mat (hidden, alphabet.size());
      dummy_em_mat.init_alphabet (alphabet);

      // now create the real matrix
      const int m = matrix_set.add_matrix (1, Rev, hidden);
      ECFG_chain& chain = matrix_set.chain[m];
      for (int c = 0; c < hidden; ++c)
	{
	  chain.class_labels[c].clear();
	  chain.class_labels[c] << c + 1;
	}
      EM_matrix_base& mat = *chain.matrix;
      mat.init_alphabet (dummy_em_mat.hidden_alphabet);

      // populate the matrix
      const double rate_mul = (double) (s+1) / (double) zones;
      const int S = alphabet.size();
      for (int c = 0; c < hidden; ++c)
	{
	  for (int i = 0; i < alphabet.size(); ++i)
	    {
	      mat.pi[c*S+i] = pam_prior[i] / (double) hidden;
	      for (int j = 0; j < alphabet.size(); ++j)
		mat.X[0](c*S+i,c*S+j) = rate_mul * pam_matrix(i,j) / (double) (c+1);
	      for (int d = 0; d < hidden; ++d)
		mat.X[0](c*S+i,d*S+i) = rate_mul / (double) S;
	    }
	}
      mat.update();

      state_info[s] = ECFG_state_info (1, 0);
      state_info[s].matrix = m;
      state_info[s].mul[0] = 1;
      const sstring state_char (1, (char) (s < 10 ? '0' + s : 'a' + s - 10));
      state_info[s].name << "S" << state_char;
      if (zones > 1)
	state_info[s].annot[sstring("State")][state_char] = PFunc(1.);
      matrix_set.chain[m].state[0].clear();
      matrix_set.chain[m].state[0] << "A" << state_char;

      transition (Start, s) = transition (Start, End) = -Prob2Score (zones + 1);
      for (int d = 0; d < zones; ++d)
	if (s != d)
	  transition (s, d) = -Prob2Score (2*zones);
      transition (s, End) = -Prob2Score (2*zones);
      transition (s, s) = -Prob2Score (2);
    }
}
