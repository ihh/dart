#include "empath/repeater.h"

Repeater::Repeater (Trainable& unit_model,
		    const int min_units,
		    const int max_units)
  : Trainable (unit_model),
    unit_model (unit_model),
    min_units (min_units),
    max_units (max_units)
{
  // create skip_unit PGroup & assign prior
  skip_unit = pscore.new_group (max_units + 1 - min_units, "No. of units to skip");
  vector<Prob> skip_pseud (max_units + 1 - min_units, 1.0);
  Dirichlet_mixture skip_dirichlet (skip_pseud);
  prior.assign (skip_unit, skip_dirichlet);
  // resize HMM
  hmm = Single_PHMM (total_states(), unit_model.hmm.alphabet());
  // make (max_units) copies of unit HMM
  for (int unit = 0; unit < max_units; ++unit)
    {
      hmm.state_type[unit_fwdpad(unit)]   = Single_state_typing::Emit;
      hmm.state_type[unit_revpad(unit)]   = Single_state_typing::Emit;
      hmm.transition (unit_fwdend(unit), unit_fwdpad(unit)) = unit_model.null_extend.YES;
      hmm.transition (unit_revend(unit), unit_revpad(unit)) = unit_model.null_extend.YES;
      hmm.transition (unit_fwdpad(unit), unit_fwdpad(unit)) = unit_model.null_extend.YES;
      hmm.transition (unit_revpad(unit), unit_revpad(unit)) = unit_model.null_extend.YES;
      for (int sym = 0; sym < hmm.alphabet().size(); ++sym)
	{
	  hmm.emit[unit_fwdpad(unit)][sym] = unit_model.null_emit[sym];
	  hmm.emit[unit_revpad(unit)][sym] = unit_model.null_emit[sym];
	}
      for (int i = 0; i < unit_model.hmm.states(); ++i)
	{
	  // set state type & emissions
	  hmm.state_type[unit_state(unit,i)] = unit_model.hmm.state_type[i];
	  hmm.emit[unit_state(unit,i)] = unit_model.hmm.emit[i];
	  // set transitions
	  for (int j = 0; j < unit_model.hmm.states(); ++j)
	    hmm.transition (unit_state(unit,i), unit_state(unit,j)) = unit_model.hmm.transition (i, j);
	  // set metascores
	  hmm.metascore_idx[unit_state(unit,i)] = unit_model.hmm.metascore_idx[i];
	}
    }
  // join consecutive units together
  for (int unit = 0; unit < max_units; ++unit)
    {
      if (unit <= max_units - min_units)
	{
	  hmm.transition (Start, unit_fwdstart(unit)) = skip_unit[unit];
	  hmm.transition (Start, unit_revstart(unit)) = skip_unit[unit];
 	}
      if (unit < max_units - 1)
	{
	  hmm.transition (unit_fwdend(unit), unit_fwdstart(unit+1)) = 1;
	  hmm.transition (unit_revend(unit), unit_revstart(unit+1)) = 1;
	  hmm.transition (unit_fwdpad(unit), unit_fwdstart(unit+1)) = 1;
	  hmm.transition (unit_revpad(unit), unit_revstart(unit+1)) = 1;
	}
    }
  if (min_units <= 0)
    hmm.transition (Start, End) = skip_unit[max_units];
  hmm.transition (unit_fwdend(max_units-1), End) = 1;
  hmm.transition (unit_revend(max_units-1), End) = 1;
  hmm.transition (unit_fwdpad(max_units-1), End) = 1;
  hmm.transition (unit_revpad(max_units-1), End) = 1;
  // output
  if (CTAGGING(3,REPEATER)) { CL << "Repeater model:\n"; hmm.show (CL); }
}

void Repeater::instruct_turtle (ostream& o) const
{
  unit_model.instruct_turtle (o);
}
