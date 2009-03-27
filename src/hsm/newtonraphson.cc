#include <algorithm>
#include "hsm/newtonraphson.h"
#include "util/math_fn.h"
#include "util/logfile.h"
#include "util/vector_output.h"

vector<double> Newton_Raphson::iterate (Grad_function& f, const vector<double>& xseed, const Constraints& xmin, double tol, int max_iter)
{
  if ((int) xseed.size() != f.dim()) THROWEXPR ("Mis-sized Newton-Raphson seed: (" << xseed << ")");
  // set up vars
  vector<double> x = xseed;
  double f_val;
  vector<double> grad_f_val (xseed.size(), (double) 0);
  double best_f = 0;
  vector<double> best_x;
  int iter = 0;
  while (1)
    {
      for_const_contents (Constraints, xmin, xm)
	x[xm->first] = max (x[xm->first], xm->second);
      f.eval (x, f_val, grad_f_val);
      if (iter == 0 || f_val < best_f)
	{
	  best_f = f_val;
	  best_x = x;
	}
      ++iter;
      if (CTAGGING(3,NEWTON_RAPHSON NEWTON_RAPHSON_GRAD))
	{
	  CL << "Newton-Raphson iteration #" << iter << ": f=" << f_val << " min=" << best_f << "\n";
	  if (CTAGGING(3,NEWTON_RAPHSON_GRAD))
	    CL << " x=(" << x << ") grad(f)=(" << grad_f_val << ")\n";
	}
      if (abs (f_val) < tol)
	{
	  CTAG(5,NEWTON_RAPHSON) << "Newton-Raphson converged after " << iter << " iterations\n";
	  break;
	}
      if (iter >= max_iter)
	{
	  CTAG(5,NEWTON_RAPHSON) << "Newton-Raphson hit " << max_iter << " iterations; using min f=" << best_f << "\n";
	  x = best_x;
	  break;
	}
      // find magnitude of gradient step
      double grad_mag_sq = 0;
      for (int i = 0; i < (int) xseed.size(); ++i)
	grad_mag_sq += grad_f_val[i] * grad_f_val[i];
      // do a step
      const double x_delta_size = f_val / grad_mag_sq;
      for (int i = 0; i < (int) xseed.size(); ++i)
	x[i] -= x_delta_size * grad_f_val[i];
    }
  return x;
}

Diff::Stack_index::Stack_index (Diff& D, bool _update) : D(D), const_i(0), param_i(0)
{
  if (_update) update();
}

void Diff::Stack_index::update()
{
  update (D.opcode.size());
}

void Diff::Stack_index::update (int new_size)
{
  const int old_size = l_op_idx.size();

  l_op_idx.reserve (new_size);
  r_op_idx.reserve (new_size);
  start_idx.reserve (new_size);
  const_idx.reserve (new_size);
  param_idx.reserve (new_size);

  for (int i = old_size; i < new_size; ++i)
    {
      int l_i = -1;
      int r_i = -1;
      int start_i = -1;
      switch (D.opcode[i])
	{
	case CONST:
	  start_i = i;
	  ++const_i;
	  break;
	case PARAM:
	  start_i = i;
	  ++param_i;
	  break;
	case ADD:
	case MUL:
	  {
	    if (idx_stack.size() < 2) THROWEXPR ("Diff contains binary operator with <2 operands");
	    r_i = idx_stack.top();
	    idx_stack.pop();
	    l_i = idx_stack.top();
	    idx_stack.pop();
	    start_i = start_idx[l_i];
	    break;
	  }
	case LOG:
	case EXP:
	  {
	    if (idx_stack.size() < 1) THROWEXPR ("Diff contains unary operator with no operand");
	    l_i = idx_stack.top();
	    idx_stack.pop();
	    start_i = start_idx[l_i];
	    break;
	  }
	case POW:
	  {
	    if (idx_stack.size() < 1) THROWEXPR ("Diff contains unary operator with no operand");
	    l_i = idx_stack.top();
	    idx_stack.pop();
	    start_i = start_idx[l_i];
	    ++const_i;
	    break;
	  }
	default:
	  THROWEXPR ("Unknown operand in Diff");
	  break;
	}
      l_op_idx.push_back (l_i);
      r_op_idx.push_back (r_i);
      start_idx.push_back (start_i);
      const_idx.push_back (const_i);
      param_idx.push_back (param_i);
      idx_stack.push (i);
    }
  if (new_size == (int) D.opcode.size() && ((int) D.constant.size() != const_i || (int) D.param.size() != param_i))
    THROW Standard_exception ("Corrupt Diff");
}

void Diff::Stack_index::assert_single_stack() const
{
  if (D.opcode.size() > 0 && idx_stack.size() != 1)
    THROWEXPR ("Diff not single-valued");
}

void Diff::Stack_index::insert_node (int expr_idx, Diff& dest, int opcode_pos, int const_pos, int param_pos) const
{
  const int expr_begin_idx = start_idx[expr_idx];
  const int const_begin_idx = expr_begin_idx < 1 ? 0 : const_idx[expr_begin_idx - 1];
  const int param_begin_idx = expr_begin_idx < 1 ? 0 : param_idx[expr_begin_idx - 1];
  
  dest.opcode.insert (dest.opcode.begin() + opcode_pos,
		      D.opcode.begin() + expr_begin_idx,
		      D.opcode.begin() + expr_idx + 1);

  dest.constant.insert (dest.constant.begin() + const_pos,
			D.constant.begin() + const_begin_idx,
			D.constant.begin() + const_idx[expr_idx]);

  dest.param.insert (dest.param.begin() + param_pos,
		     D.param.begin() + param_begin_idx,
		     D.param.begin() + param_idx[expr_idx]);
}

void Diff::Stack_index::append_node (int expr_idx, Diff& dest) const
{
  insert_node (expr_idx, dest, dest.opcode.size(), dest.constant.size(), dest.param.size());
}

void Diff::Stack_index::delete_node (int expr_idx)
{
  const int expr_begin_idx = start_idx[expr_idx];
  const int const_begin_idx = expr_begin_idx < 1 ? 0 : const_idx[expr_begin_idx - 1];
  const int param_begin_idx = expr_begin_idx < 1 ? 0 : param_idx[expr_begin_idx - 1];

  D.opcode.erase (D.opcode.begin() + expr_begin_idx, D.opcode.begin() + expr_idx + 1);
  D.constant.erase (D.constant.begin() + const_begin_idx, D.constant.begin() + const_idx[expr_idx]);
  D.param.erase (D.param.begin() + param_begin_idx, D.param.begin() + param_idx[expr_idx]);

  forget (expr_begin_idx);
}

void Diff::Stack_index::swap_last_two_nodes()
{
  idx_stack.pop();
  const int l = idx_stack.top();
  idx_stack.pop();
  const int opcode_i = l_op_idx.size();
  insert_node (l, D, opcode_i, const_i, param_i);
  delete_node (l);
  update (opcode_i);
}

void Diff::Stack_index::forget (int pos)
{
  l_op_idx.erase (l_op_idx.begin() + pos, l_op_idx.end());
  r_op_idx.erase (r_op_idx.begin() + pos, r_op_idx.end());
  start_idx.erase (start_idx.begin() + pos, start_idx.end());
  const_idx.erase (const_idx.begin() + pos, const_idx.end());
  param_idx.erase (param_idx.begin() + pos, param_idx.end());

  const_i = pos > 0 ? const_idx.back() : 0;
  param_i = pos > 0 ? param_idx.back() : 0;
}

sstring Diff::make_param_name (int i) const
{
  if (param_name) return (*param_name)[i];
  sstring n;
  n << "X" << i;
  return n;
}

int Diff::dim() const
{
  if (param.size() == 0) return 0;
  return *max_element (param.begin(), param.end()) + 1;
}

set<int> Diff::params() const
{
  return set<int> (param.begin(), param.end());
}

double Diff::eval (const vector<double>& x) const
{
  if (is_null()) return 0;
  vector<double>::const_iterator const_ptr = constant.begin();
  vector<int>::const_iterator param_ptr = param.begin();
  vector<double> val;
  val.reserve (opcode.size());
  for (int i = 0; i < (int) opcode.size(); ++i)
    switch (opcode[i])
      {
      case CONST:
	val.push_back (*const_ptr++);
	break;
      case PARAM:
	val.push_back (x[*param_ptr++]);
	break;
      case ADD:
	{
	  const double r_val = val.back();
	  val.pop_back();
	  val.back() += r_val;
	  break;
	}
      case MUL:
	{
	  const double r_val = val.back();
	  val.pop_back();
	  val.back() *= r_val;
	  break;
	}
      case LOG:
	val.back() = Math_fn::math_log (val.back());
	break;
      case EXP:
	val.back() = Math_fn::math_exp (val.back());
	break;
      case POW:
	val.back() = Math_fn::math_pow (val.back(), *const_ptr++);
	break;
      default:
	THROWEXPR ("When is Spain?");
	break;
      }
  if (val.size() != 1) THROWEXPR ("Bastardo! Diff stack empty");
  return val.back();
}

Diff Diff::derivative (int deriv_idx) const
{
  Stack_index s ((Diff&) *this);
  s.assert_single_stack();
  Diff deriv;
  deriv.param_name = param_name;
  Stack_index sd (deriv);
  vector<int> deriv_expr;  // mapping from nodes to their derivatives
  for (int i = 0; i < (int) opcode.size(); ++i)
    {
      const int l_op = s.l_op_idx[i];
      const int r_op = s.r_op_idx[i];
      const int l_op_deriv = l_op < 0 ? l_op : deriv_expr[l_op];
      const int r_op_deriv = r_op < 0 ? r_op : deriv_expr[r_op];
      bool vanished = 0;
      const bool l_vanished = l_op_deriv < 0;
      const bool r_vanished = r_op_deriv < 0;
      switch (opcode[i])
	{
	case CONST:
	  vanished = 1;
	  break;
	case PARAM:
	  if (param [s.param_idx[i] - 1] != deriv_idx)
	    vanished = 1;
	  else
	    {
	      deriv.constant.push_back (1.0);
	      deriv.opcode.push_back (CONST);
	    }
	  break;
	case ADD:
	  if (l_vanished && r_vanished)
	    vanished = 1;
	  else if (!l_vanished && !r_vanished)
	    deriv.opcode.push_back (ADD);
	  break;
	case MUL:
	  if (l_vanished && r_vanished)
	    vanished = 1;
	  else
	    {
	      // need to go from "f' g'" to "g' f MUL f' g MUL ADD"
	      if (!r_vanished)
		{
		  // "g'" ...
		  s.append_node (l_op, deriv);  // -> "g' f"
		  deriv.opcode.push_back (MUL);  // -> "g' f MUL"
		  sd.update();
		}
	      if (!l_vanished)
		{
		  // "f' g' f MUL" ...
		  if (!r_vanished) sd.swap_last_two_nodes();  // -> "g' f MUL f'"
		  s.append_node (r_op, deriv);  // -> "g' f MUL f' g"
		  deriv.opcode.push_back (MUL);  // -> "g' f MUL f' g MUL"
		  if (!r_vanished) deriv.opcode.push_back (ADD);  // -> "g' f MUL f' g MUL ADD"
		}
	    }
	  break;
	case LOG:
	  if (l_vanished)
	    vanished = 1;
	  else
	    {
	      // need to go from "f'" to "f' f POW(-1) MUL"
	      s.append_node (l_op, deriv);
	      deriv.opcode.push_back (POW);
	      deriv.constant.push_back (-1);
	      deriv.opcode.push_back (MUL);
	    }
	  break;
	case EXP:
	  if (l_vanished)
	    vanished = 1;
	  else
	    {
	      // need to go from "f'" to "f' f EXP MUL"
	      s.append_node (l_op, deriv);
	      deriv.opcode.push_back (EXP);
	      deriv.opcode.push_back (MUL);
	    }
	  break;
	case POW:
	  if (l_vanished)
	    vanished = 1;
	  else
	    {
	      const double c = constant [s.const_idx[i] - 1];
	      if (c != 1)
		{
		  // need to go from "f'" to "f' f POW(c-1) CONST(c) MUL MUL"
		  deriv.opcode.push_back (CONST);
		  deriv.constant.push_back (c);
		  s.append_node (l_op, deriv);
		  if (c != 2)
		    {
		      deriv.opcode.push_back (POW);
		      deriv.constant.push_back (c - 1);
		    }
		  deriv.opcode.push_back (MUL);
		  deriv.opcode.push_back (MUL);
		}
	    }
	  break;
	default:
	  THROWEXPR ("Why is Hitler?");
	  break;
	}
      if (vanished)
	deriv_expr.push_back (-1);
      else
	deriv_expr.push_back (deriv.opcode.size() - 1);
      sd.update();
    }
  deriv.simplify();
  return deriv;
}

sstring Diff::to_string() const
{
  if (is_null()) return sstring ("0");
  Stack_index stack_idx ((Diff&) *this);
  stack<sstring> text_stack;
  for (int i = 0; i < (int) opcode.size(); ++i)
    {
      sstring text;
      const Opcode l_op = stack_idx.l_op_idx[i] < 0 ? NOP : opcode[stack_idx.l_op_idx[i]];
      const Opcode r_op = stack_idx.r_op_idx[i] < 0 ? NOP : opcode[stack_idx.r_op_idx[i]];
      switch (opcode[i])
	{
	case CONST:
	  text << constant[stack_idx.const_idx[i] - 1];
	  break;

	case PARAM:
	  if (param_name != 0)
	    text << (*param_name)[param[stack_idx.param_idx[i] - 1]];
	  else
	    text << 'X' << param[stack_idx.param_idx[i] - 1];
	  break;

	case ADD:
	  {
	    const sstring r_text = text_stack.top();
	    text_stack.pop();
	    const sstring l_text = text_stack.top();
	    text_stack.pop();
	    text << l_text;
	    if (r_text[0] != '-')
	      text << '+';
	    text << r_text;
	    break;
	  }

	case MUL:
	  {
	    const sstring r_text = text_stack.top();
	    text_stack.pop();
	    const sstring l_text = text_stack.top();
	    text_stack.pop();
	    if (l_op == ADD)
	      text << '(' << l_text << ')';
	    else
	      text << l_text;
	    text << '*';
	    if (r_op == ADD || r_text[0] == '-')
	      text << '(' << r_text << ')';
	    else
	      text << r_text;
	    break;
	  }

	case LOG:
	  {
	    const sstring l_text = text_stack.top();
	    text_stack.pop();
	    text << "log(" << l_text << ')';
	    break;
	  }

	case EXP:
	  {
	    const sstring l_text = text_stack.top();
	    text_stack.pop();
	    text << "exp(" << l_text << ')';
	    break;
	  }

	case POW:
	  {
	    const sstring l_text = text_stack.top();
	    text_stack.pop();
	    if (l_op == MUL || l_op == ADD)
	      text << '(' << l_text << ')';
	    else
	      text << l_text;
	    text << '^';
	    const double c = constant[stack_idx.const_idx[i] - 1];
	    if (c < 0)
	      text << '(' << c << ')';
	    else text << c;
	    break;
	  }

	default:   // never get here
	  THROWEXPR ("Where are the Snowdens of yesteryear?");
	  break;
	}
      text_stack.push (text);
    }
  const sstring result = text_stack.top();
  return result;
}

sstring Diff::params_to_string (const vector<double>& x) const
{
  const int N = dim();
  sstring ps;
  for (int i = 0; i < N; ++i)
    {
      if (param_name) ps << (*param_name)[i];
      else ps << 'X' << i;
      ps << '=' << x[i];
      if (i < N - 1) ps << ' ';
    }
  return ps;
}

void Diff::append (const Diff& diff)
{
  opcode.insert (opcode.end(), diff.opcode.begin(), diff.opcode.end());
  constant.insert (constant.end(), diff.constant.begin(), diff.constant.end());
  param.insert (param.end(), diff.param.begin(), diff.param.end());
}

bool Diff::is_null() const
{
  return opcode.size() == 0;
}

void Diff::clear()
{
  opcode.clear();
  constant.clear();
  param.clear();
}

Diff& Diff::operator*= (const Diff& diff)
{
  if (!is_null() && !diff.is_null())
  {
    append (diff);
    opcode.push_back (MUL);
  }
  return *this;
}

Diff& Diff::operator+= (const Diff& diff)
{
  const bool lhs_null = is_null();
  if (!diff.is_null())
  {
    append (diff);
    if (!lhs_null) opcode.push_back (ADD);
  }
  return *this;
}

Diff& Diff::operator-= (const Diff& diff)
{
  const bool lhs_null = is_null();
  if (!diff.is_null())
  {
    opcode.push_back (CONST);
    constant.push_back (-1);
    append (diff);
    opcode.push_back (MUL);
    if (!lhs_null) opcode.push_back (ADD);
  }
  return *this;
}

Diff& Diff::operator/= (const Diff& diff)
{
  if (diff.is_null()) THROWEXPR ("Symbolic division by zero");
  if (!is_null())
  {
    append (diff);
    opcode.push_back (POW);
    constant.push_back (-1);
    opcode.push_back (MUL);
  }
  return *this;
}

Diff operator* (const Diff& a, const Diff& b)
{
  Diff diff (a);
  diff *= b;
  return diff;
}

Diff operator+ (const Diff& a, const Diff& b)
{
  Diff diff (a);
  diff += b;
  return diff;
}

Diff operator- (const Diff& a, const Diff& b)
{
  Diff diff (a);
  diff -= b;
  return diff;
}

Diff operator/ (const Diff& a, const Diff& b)
{
  Diff diff (a);
  diff /= b;
  return diff;
}

Diff Diff::log() const
{
  Diff diff (*this);
  diff.opcode.push_back (LOG);
  return diff;
}

Diff Diff::exp() const
{
  Diff diff (*this);
  diff.opcode.push_back (EXP);
  return diff;
}

Diff Diff::pow (double p) const
{
  Diff diff (*this);
  diff.opcode.push_back (POW);
  diff.constant.push_back (p);
  return diff;
}

Diff Diff::sum (vector<Diff>& series)
{
  Diff diff;
  for_const_contents (vector<Diff>, series, s) diff += *s;
  diff.simplify();
  return diff;
}

void Diff::simplify()
{
  if (CTAGGING(3,DIFF_SIMPLIFY)) CL << "Simplifying " << to_string() << "\n";
  // find MUL/ADD cliques; replace multiple leaves (CONST's) with NOP's, propagating upwards
  Stack_index stack_idx (*this);
  vector<int> parent (opcode.size(), (int) -1);
  vector<int> clique_root (opcode.size());
  vector<int> clique_leaf (opcode.size(), (int) -1);
  for (int i = opcode.size() - 1; i >= 0; --i)
    {
      const Opcode op = opcode[i];
      int& root = clique_root[i];

      if (parent[i] < 0) root = i;
      else root = clique_root[parent[i]];  // may change if op is ADD or MUL; see below

      if (op == ADD || op == MUL)
	{
	  if (opcode[root] != op) root = i;  // keep ADD & MUL cliques distinct
	  parent[stack_idx.l_op_idx[i]] = i;
	  parent[stack_idx.r_op_idx[i]] = i;
	}
      else if (op == CONST)
	{
	  int& leaf = clique_leaf[root];
	  if (leaf < 0)
	    leaf = i;  // first leaf of clique
	  else
	    {
	      // fold leaves
	      const double node_const = constant[stack_idx.const_idx[i] - 1];
	      double& leaf_const = constant[stack_idx.const_idx[leaf] - 1];
	      const Opcode root_op = opcode[root];
	      if (root_op == ADD)  // ADD
		leaf_const += node_const;
	      else  // MUL
		leaf_const *= node_const;
	      // replace CONST with NOP
	      int nop_node = i;
	      opcode[nop_node] = NOPC;
	      // propagate NOP upwards
	      while (opcode[nop_node] != root_op) nop_node = parent[nop_node];
	      opcode[nop_node] = NOP;
	    }
	}
    }
  // input
  vector<Opcode>::const_iterator opcode_i = opcode.begin();
  vector<double>::const_iterator const_i = constant.begin();
  vector<int>::const_iterator param_i = param.begin();
  // output
  vector<Opcode>::iterator opcode_o = opcode.begin();
  vector<double>::iterator const_o = constant.begin();
  vector<int>::iterator param_o = param.begin();
  // merge "POW" runs; excise "LOG EXP", "EXP LOG", "NOP"
  while (opcode_i < opcode.end())
    {
      const Opcode op = *opcode_i++;
      switch (op)
	{
	case CONST:
	  *opcode_o++ = op;
	  *const_o++ = *const_i++;
	  break;
	case PARAM:
	  *opcode_o++ = op;
	  *param_o++ = *param_i++;
	  break;
	case POW:
	  {
	    double p = *const_i++;
	    while (opcode_i < opcode.end())
	      {
		if (*opcode_i != POW) break;
		p *= *const_i++;
		++opcode_i;
	      }
	    if (p != 0)
	      {
		*opcode_o++ = op;
		*const_o++ = p;
	      }
	    break;
	  }
	case LOG:
	  if (opcode_i < opcode.end() ? *opcode_i == EXP : 0)
	    ++opcode_i;
	  else
	    *opcode_o++ = op;
	  break;
	case EXP:
	  if (opcode_i < opcode.end() ? *opcode_i == LOG : 0)
	    ++opcode_i;
	  else
	    *opcode_o++ = op;
	  break;
	case ADD:
	case MUL:
	  *opcode_o++ = op;
	  break;
	case NOPC:
	  ++const_i;
	  break;
	case NOPP:
	  ++param_i;
	  break;
	case NOP:
	  break;
	default:
	  THROWEXPR ("WTF?");
	  break;
	}
    }
  const bool simplified = opcode_o < opcode_i;
  opcode.erase (opcode_o, opcode.end());
  constant.erase (const_o, constant.end());
  param.erase (param_o, param.end());
  if (simplified && CTAGGING(3,DIFF_SIMPLIFY)) CL << "Simplified to " << to_string() << "\n";
}

Diff_const::Diff_const (double c) : Diff()
{
  opcode.push_back (CONST);
  constant.push_back (c);
}

Diff_param::Diff_param (int p) : Diff()
{
  opcode.push_back (PARAM);
  param.push_back (p);
}

Diff_log::Diff_log (const Diff& diff) : Diff (diff.log()) { }
Diff_exp::Diff_exp (const Diff& diff) : Diff (diff.exp()) { }
Diff_power::Diff_power (const Diff& diff, double p) : Diff (diff.pow (p)) { }

Laplacian::Laplacian (const Diff& D) : D(D), d2D_dxdy (D.dim(), Diff_map())
{
  if (CTAGGING(3,LAPLACE LAPLACE_DERIV)) CL << "Finding grad^2(D) where D = " << D.to_string() << "\n";
  for (int i = 0; i < D.dim(); ++i)
    {
      dD_dx.push_back (D.derivative (i));
      if (CTAGGING(3,LAPLACE_DERIV))
	CL << "diff(D," << D.make_param_name (i) << ") = " << dD_dx[i].to_string() << "\n";
      for_tmp_contents (set<int>, dD_dx[i].params(), j)
	if (*j >= i)
	  {
	    d2D_dxdy[i].insert (Diff_map::value_type (*j, dD_dx[i].derivative (*j)));
	    if (CTAGGING(3,LAPLACE_DERIV))
	      {
		CL << "diff(D,";
		if (D.param_name) CL << (*D.param_name)[i]; else CL << "X" << i;
		CL << ",";
		if (D.param_name) CL << (*D.param_name)[*j]; else CL << "X" << *j;
		CL << ") = " << d2D_dxdy[i][*j].to_string() << "\n";
	      }
	  }
    }
}

int Laplacian::dim() { return D.dim(); }

void Laplacian::eval (const vector<double>& x, double& f, vector<double>& grad_f)
{
  const int dim_D = D.dim();
  vector<double> d (dim_D);
  f = 0;
  grad_f = vector<double> (dim_D, (double) 0);
  int i_max = -1;
  double d_max = 0;
  for (int i = 0; i < dim_D; ++i)
    {
      const double d_i = dD_dx[i].eval (x);
      const double d_sq = d_i * d_i;
      if (i == 0 || d_sq > d_max) { d_max = d_sq; i_max = i; }
      d[i] = d_i;
      f += d_sq;
    }
  if (CTAGGING(3,LAPLACE LAPLACE_PARAMS))
    {
      CL << "D = " << D.eval (x) << ", diff(D," << D.make_param_name (i_max) << ") dominates\n";
      if (CTAGGING(2,LAPLACE_PARAMS))
	CL << " (" << D.params_to_string (x) << ")\n";
    }
  for (int i = 0; i < dim_D; ++i)
    for_const_contents (Diff_map, d2D_dxdy[i], j_d2D)
      {
	const int j = j_d2D->first;
	const double d2_ij = j_d2D->second.eval (x);
	grad_f[j] += 2 * d[i] * d2_ij;
	if (i != j)
	  grad_f[i] += 2 * d[j] * d2_ij;
      }
}
