#include <algorithm>
#include "seq/pfunc.h"

bool PFunc::is_null() const
{
  return opcode.size() == 0;
}

bool PFunc::is_zero() const
{
  if (!is_null())
    if (opcode.size() == 1)
      if (opcode[0] == CONST)
	if (const_sc[0] == -InfinityScore)
	  return TRUE;
  return FALSE;
}

bool PFunc::is_one() const
{
  if (!is_null())
    if (opcode.size() == 1)
      if (opcode[0] == CONST)
	if (const_sc[0] == 0)
	  return TRUE;
  return FALSE;
}

set<PVar> PFunc::vars() const
{
  return set<PVar> (pvar.begin(), pvar.end());
}

set<PGroup> PFunc::groups() const
{
  set<PGroup> g;
  for_const_contents (vector<PVar>, pvar, p) g.insert ((*p).group_idx);
  return g;
}

PFunc::Stack_index::Stack_index (const PFunc& f)
  : f (f),
    l_operand_idx (f.opcode.size()),
    r_operand_idx (f.opcode.size()),
    cp_idx (f.opcode.size()),
    node_sc (f.opcode.size())
{
  // test for null funcs
  if (f.is_null())
    THROW Standard_exception ("Tried to evaluate a null PFunc");

  // set up the pointers
  int const_idx = 0;
  int var_idx = 0;
  vector<int> idx_stack;
  for (int i = 0; i < (int) f.opcode.size(); ++i)
    {
      switch (f.opcode[i])
	{
	case CONST:
	  cp_idx[i] = const_idx++;
	  idx_stack.push_back(i);
	  break;
	case PVAR:
	case INERT_PVAR:
	  cp_idx[i] = var_idx++;
	  idx_stack.push_back(i);
	  break;
	case ADD:
	case MUL:
	case DIV:
	  if (idx_stack.size() < 2) THROW Standard_exception ("PFunc contains binary operator with <2 operands");
	  r_operand_idx[i] = idx_stack.back();
	  idx_stack.pop_back();
	  l_operand_idx[i] = idx_stack.back();
	  idx_stack.pop_back();
	  idx_stack.push_back(i);
	  break;
	case POW:
	  if (idx_stack.size() < 1) THROW Standard_exception ("PFunc contains POW operator with no function to exponentiate");
	  cp_idx[i] = const_idx++;
	  l_operand_idx[i] = idx_stack.back();
	  idx_stack.pop_back();
	  idx_stack.push_back(i);
	  break;
	default:
	  THROW Standard_exception ("Unknown operand in PFunc");
	  break;
	}
    }
  if ((int) idx_stack.size() != 1 || (int) f.const_sc.size() != const_idx || (int) f.pvar.size() != var_idx)
    THROW Standard_exception ("Corrupt PFunc");
}

void PFunc::Stack_index::eval_node_sc (const PScores& var_sc)
{
  for (int i = 0; i < (int) f.opcode.size(); ++i)
    {
      switch (f.opcode[i])
	{
	case CONST:
	  node_sc[i] = f.const_sc[cp_idx[i]];
	  break;
	case POW:
	  node_sc[i] = (Score) (FScore2Bits (f.const_sc[cp_idx[i]]) * (double) node_sc[l_operand_idx[i]]);
	  break;
	case PVAR:
	case INERT_PVAR:
	  node_sc[i] = var_sc [(f.pvar[cp_idx[i]])];
	  break;
	case ADD:
	  node_sc[i] = ScorePSum (node_sc[l_operand_idx[i]], node_sc[r_operand_idx[i]]);
	  break;
	case MUL:
	  node_sc[i] = ScorePMul (node_sc[l_operand_idx[i]], node_sc[r_operand_idx[i]]);
	  break;
	case DIV:
	  node_sc[i] = ScorePMul (node_sc[l_operand_idx[i]], -node_sc[r_operand_idx[i]]);
	  break;
	default:   // never get here
	  break;
	}
    }
}

Score PFunc::eval_sc (const PScores& var_sc) const
{
  Stack_index stack_idx (*this);
  stack_idx.eval_node_sc (var_sc);
  return stack_idx.node_sc.back();
}

Prob PFunc::eval_dlogf_dlogp (const PScores& var_sc, const PVar& pv) const
{
  Stack_index stack_idx (*this);
  stack_idx.eval_node_sc (var_sc);
  vector<Prob> dlognodesc_dlogp (opcode.size());
  for (int i = 0; i < (int) opcode.size(); ++i)
    {
      switch (opcode[i])
	{
	case CONST:
	case INERT_PVAR:
	  dlognodesc_dlogp[i] = 0.0;
	  break;
	case PVAR:
	  dlognodesc_dlogp[i] = pvar[stack_idx.cp_idx[i]] == pv ? 1.0 : 0.0;
	  break;
	case ADD:
	  {
	    const int l_idx = stack_idx.l_operand_idx[i];
	    const int r_idx = stack_idx.r_operand_idx[i];
	    const Prob lop = Score2Prob(stack_idx.node_sc[l_idx]);
	    const Prob rop = Score2Prob(stack_idx.node_sc[r_idx]);
	    const Prob dlop_dlogp = lop * dlognodesc_dlogp[l_idx];
	    const Prob drop_dlogp = rop * dlognodesc_dlogp[r_idx];
	    dlognodesc_dlogp[i] = (dlop_dlogp + drop_dlogp) / (lop + rop);
	  }
	  break;
	case MUL:
	  dlognodesc_dlogp[i] = dlognodesc_dlogp[stack_idx.l_operand_idx[i]] + dlognodesc_dlogp[stack_idx.r_operand_idx[i]];
	  break;
	case DIV:
	  dlognodesc_dlogp[i] = dlognodesc_dlogp[stack_idx.l_operand_idx[i]] - dlognodesc_dlogp[stack_idx.r_operand_idx[i]];
	  break;
	case POW:
	  dlognodesc_dlogp[i] = (Score) (FScore2Bits (const_sc[stack_idx.cp_idx[i]]) * (double) dlognodesc_dlogp[stack_idx.l_operand_idx[i]]);
	  break;
	default:   // never get here
	  break;
	}
    }
  return dlognodesc_dlogp.back();
}

void PFunc::inc_var_counts (PCounts& var_counts, const PScores& var_sc, Prob func_count, Prob func_wait) const
{
  set<PVar> pvars = vars();
  const Prob f = Score2Prob (eval_sc (var_sc));
  for_const_contents (set<PVar>, pvars, pvar)
    {
      const Prob dlogf_dlogp = eval_dlogf_dlogp (var_sc, *pvar);
      const Prob p = Score2Prob (var_sc[*pvar]);
      const Prob df_dp = dlogf_dlogp * f / p;

      var_counts [(*pvar)] += dlogf_dlogp * func_count;
      if (var_counts.group_size (pvar->group_idx) == 1)  // rate variable?
	var_counts.wait[pvar->group_idx] += df_dp * func_wait;
    }
}

PFunc::PFunc() : opcode(), const_sc(), pvar() { }

PFunc::PFunc (Prob constant)
  : opcode (1, CONST),
    const_sc (1, (FScore) Prob2FScore(constant)),
    pvar()
{ }

PFunc::PFunc (PFunc f, Prob exponent)
  : opcode (f.opcode),
    const_sc (f.const_sc),
    pvar (f.pvar)
    
{
  if (f.is_null())
    THROW Standard_exception ("Tried to exponentiate a null PFunc");

  opcode.push_back (POW);
  const_sc.push_back ((FScore) Bits2FScore(exponent));
}

PFunc::PFunc (PVar pv) : opcode (1, (Opcode) PVAR), const_sc(), pvar (1, pv) { }

PFunc::PFunc (Opcode op, PFunc lfunc, PFunc rfunc)
  : opcode (lfunc.opcode.size() + rfunc.opcode.size() + 1),
    const_sc (lfunc.const_sc.size() + rfunc.const_sc.size()),
    pvar (lfunc.pvar.size() + rfunc.pvar.size())
{
  if (lfunc.is_null() || rfunc.is_null())
    THROW Standard_exception ("Tried to use a null PFunc in a binary context");

  copy (lfunc.opcode.begin(),   lfunc.opcode.end(),   opcode.begin());
  copy (lfunc.const_sc.begin(), lfunc.const_sc.end(), const_sc.begin());
  copy (lfunc.pvar.begin(),     lfunc.pvar.end(),     pvar.begin());

  copy (rfunc.opcode.begin(),   rfunc.opcode.end(),   opcode.begin()   + lfunc.opcode.size());
  copy (rfunc.const_sc.begin(), rfunc.const_sc.end(), const_sc.begin() + lfunc.const_sc.size());
  copy (rfunc.pvar.begin(),     rfunc.pvar.end(),     pvar.begin()     + lfunc.pvar.size());

  opcode.back() = op;
}

void PFunc::make_inert()
{
  for_contents (vector<Opcode>, opcode, op)
    if (*op == PVAR)
      *op = INERT_PVAR;
}

void PFunc::show (ostream& o, const vector<vector<sstring> >* group_suffix, bool no_group_prefix) const
{
  if (is_null())
    {
      o << "[null]";
      return;
    }
  Stack_index stack_idx (*this);
  vector<sstring> node_text (opcode.size());
  for (int i = 0; i < (int) opcode.size(); ++i)
    {
      sstring& text = node_text[i];
      const int l_idx = stack_idx.l_operand_idx[i];
      const int r_idx = stack_idx.r_operand_idx[i];
      const int c_idx = stack_idx.cp_idx[i];
      switch (opcode[i])
	{
	case CONST:
	  text << FScore2Prob (const_sc[c_idx]);
	  break;

	case PVAR:
	  pvar[stack_idx.cp_idx[i]].show (text, group_suffix, no_group_prefix);
	  break;

	case INERT_PVAR:
	  text << "(# ";
	  pvar[stack_idx.cp_idx[i]].show (text, group_suffix, no_group_prefix);
	  text << ')';
	  break;

	case ADD:
	  text << node_text[l_idx] << " + " << node_text[r_idx];
	  break;

	case MUL:
	  if (opcode[l_idx] == ADD)
	    text << '(' << node_text[l_idx] << ')';
	  else
	    text << node_text[l_idx];
	  text << " * ";
	  if (opcode[r_idx] == ADD)
	    text << '(' << node_text[r_idx] << ')';
	  else
	    text << node_text[r_idx];
	  break;

	case DIV:
	  if (opcode[l_idx] == ADD)
	    text << '(' << node_text[l_idx] << ')';
	  else
	    text << node_text[l_idx];
	  text << " / ";
	  if (opcode[r_idx] == ADD || opcode[r_idx] == MUL || opcode[r_idx] == DIV)
	    text << '(' << node_text[r_idx] << ')';
	  else
	    text << node_text[r_idx];
	  break;

	case POW:
	  text << "((" << node_text[l_idx] << ") ^ " << FScore2Bits (const_sc[c_idx]) << ')';
	  break;

	default:   // never get here
	  break;
	}
    }
  o << node_text.back();
}
