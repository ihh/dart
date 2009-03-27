#include "empath/p53.h"
#include "util/math_fn.h"

P53_model::P53_model (int unit_sz, int max_spacer, PScores& pscore) :
  Trainable (4 * unit_sz + max_spacer, DNA_alphabet, pscore),
  unit_sz (unit_sz),
  max_spacer (max_spacer)
{
  // declare the PGroups
  for (int i = 0; i < unit_sz; ++i)
    {
      sstring gname;
      gname << "Match emit #" << i+1;
      match_emit.push_back (new_emit_group (gname.c_str()));
    }
  space_len = pscore.new_group (max_spacer + 1, "Space between palindromes");
  // set up the emissions
  for (int pos = 0; pos < unit_sz; ++pos)
    for (int sym = 0; sym < hmm.alphabet().size(); ++sym)
      {
	const int comp = hmm.alphabet().complement (sym);
	hmm.emit[match_state (0, pos)][sym] = hmm.emit[match_state (2, pos)][sym] = match_emit[pos][sym];
	hmm.emit[match_state (1, pos)][sym] = hmm.emit[match_state (3, pos)][sym] = match_emit[pos][comp];
      }
  for (int space_pos = 0; space_pos < max_spacer; ++space_pos)
    for (int sym = 0; sym < hmm.alphabet().size(); ++sym)
      hmm.emit[spacer_state (space_pos)][sym] = null_emit[sym];
  // set up the transitions
  // first, transitions between consecutive states in same conceptual block
  for (int pos = 0; pos < unit_sz - 1; ++pos)
    {
      hmm.transition (match_state (0, pos),     match_state (0, pos + 1)) = 1;
      hmm.transition (match_state (1, pos + 1), match_state (1, pos))     = 1;
      hmm.transition (match_state (2, pos),     match_state (2, pos + 1)) = 1;
      hmm.transition (match_state (3, pos + 1), match_state (3, pos))     = 1;
    }
  for (int space_pos = 0; space_pos < max_spacer - 1; ++space_pos)
    hmm.transition (spacer_state (space_pos), spacer_state (space_pos + 1)) = 1;
  // next, easy transitions between blocks
  hmm.transition (Start,                      match_state(0,0))         = 1;
  hmm.transition (match_state(0,unit_sz-1),   match_state(1,unit_sz-1)) = 1;
  hmm.transition (spacer_state(max_spacer-1), match_state(2,0))         = 1;
  hmm.transition (match_state(2,unit_sz-1),   match_state(3,unit_sz-1)) = 1;
  hmm.transition (match_state(3,0),           End)                      = 1;
  // finally, transitions into spacer block, parameterised by 'space' PGroup
  hmm.transition (match_state(1,0), match_state(2,0)) = space_len[0];
  for (int n_spaces = 1; n_spaces <= max_spacer; ++n_spaces)
    hmm.transition (match_state(1,0), spacer_state (max_spacer - n_spaces)) = space_len[n_spaces];
  // set up the prior for transitions into spacer block (emit prior is handled by Trainable)
  Dirichlet_mixture space_mixture (vector<Prob> (space_len.group_size, 1.0));
  prior.assign (space_len, space_mixture);
  // print ourselves out to the logfiles
  if (CTAGGING(2,P53_DISPLAY)) hmm.show(CL);
}
