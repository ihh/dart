#ifndef REPEATER_INCLUDED
#define REPEATER_INCLUDED

#include "empath/trainer.h"

class Repeater : public Trainable
{
private:
  Trainable& unit_model;
public:
  int min_units;
  int max_units;
  int min_spacers;
  int max_spacers;
  PGroup skip_unit;
  // constructor
  Repeater (Trainable& unit_model,
	    const int min_units,
	    const int max_units);
  // override virtual turtle method
  void instruct_turtle (ostream& o) const;
  // methods to get state indices
  int unit_offset   (int unit) const { return unit * (unit_model.hmm.states() + 2); }
  int unit_fwdstart (int unit) const { return unit_state (unit, Trainable::FwdStart); }
  int unit_revstart (int unit) const { return unit_state (unit, Trainable::RevStart); }
  int unit_fwdend   (int unit) const { return unit_state (unit, Trainable::FwdEnd); }
  int unit_revend   (int unit) const { return unit_state (unit, Trainable::RevEnd); }
  int unit_state    (int unit, int state) const { return unit_offset (unit) + state; }
  int unit_fwdpad   (int unit) const { return unit_offset (unit + 1) - 2; }
  int unit_revpad   (int unit) const { return unit_offset (unit + 1) - 1; }
  int total_states() const { return unit_offset (max_units); }
};

#endif /* REPEATER_INCLUDED */
