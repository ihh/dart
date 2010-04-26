#include <glob.h>

#include "util/svisitor.h"
#include "util/logfile.h"
#include "util/vector_output.h"

SExpr_macro_aliases::SExpr_macro_aliases()
{

#define SHORT(X) short2long[SEXPR_SHORTHAND_ ## X] = SEXPR_ ## X ;

  SHORT(CONCATENATE);
  SHORT(SUM);
  SHORT(MULTIPLY);
  SHORT(DIVIDE);
  SHORT(SUBTRACT);
  SHORT(MODULUS);
  SHORT(CONDITIONAL);
  SHORT(EQUALS);
  SHORT(NOT_EQUALS);
  SHORT(GREATER);
  SHORT(LESS);
  SHORT(GEQ);
  SHORT(LEQ);

}

sstring SExpr_macro_aliases::expand (const sstring& op)
{
  sstring expanded;
  const map<sstring,sstring>::const_iterator alias = short2long.find(op);
  if (alias == short2long.end())
    expanded = op;
  else
    {
      expanded = alias->second;
      if (warned.find(op) == warned.end())
	{
	  CLOGERR << "Preprocessor warning: Macro shorthand (" << op << " ARGS) can be confusing; consider (" << expanded << " ARGS) instead\n";
	  warned[op] = 1;
	}
    }
  return expanded;
}

SExpr_macro_aliases sexpr_macro_aliases;

void SExpr_visitor::preorder_visit (SExpr& sexpr)
{
  log_visit (sexpr);
  for_contents (list<SExpr>, sexpr.child, c)
    preorder_visit (*c);
}

void SExpr_visitor::postorder_visit (SExpr& sexpr)
{
  for_contents (list<SExpr>, sexpr.child, c)
    postorder_visit (*c);
  log_visit (sexpr);
}

void SExpr_visitor::log_visit (SExpr& sexpr)
{
  // log pre-visit
  if (CTAGGING(-1,SEXPR_MACROS))
    CL << "Visiting " << sexpr << '\n';

  // do the visit
  visit (sexpr);

  // log post-visit
  if (CTAGGING(-2,SEXPR_MACROS))
    CL << "Post-visit: " << sexpr << '\n';
}

void SExpr_file_operations::visit (SExpr& parent_sexpr)
{
  // prepare list of include statement nodes to erase
  list<SExprIter> erase_pos;

  // follow includes
  for_iterator (SExprIter, child_iter, parent_sexpr.child.begin(), parent_sexpr.child.end())
    {
      SExpr& sexpr = *child_iter;
      if (sexpr.is_list() && !sexpr.child.empty() && sexpr[0].is_atom())
	{
	  sstring op (sexpr[0].atom);  // make a copy, not a reference, since we will be tampering with the expression
	  CTAG(2,SEXPR_MACROS) << "Considering " << op << '\n';
	  op = sexpr_macro_aliases.expand(op);
	  // functions with one atomic argument
	  if (op == SEXPR_INCLUDE)
	    {
	      const vector<sstring> raw_includes = sexpr.atoms_to_strings (1);
	      vector<sstring> globbed_includes;

	      // expand globs
	      for_const_contents (vector<sstring>, raw_includes, raw_include)
		{
		  glob_t glob_str;
		  glob (raw_include->c_str(), GLOB_BRACE | GLOB_TILDE, (int (*)(const char *, int)) 0, &glob_str);
		  for (unsigned int g = 0; g < glob_str.gl_pathc; ++g)
		    {
		      const char* const glob_path = glob_str.gl_pathv[g];
		      CTAG(1,SEXPR_MACROS) << " Path #" << g+1 << " from glob expansion of " << *raw_include << " is " << glob_path << '\n';
		      globbed_includes.push_back (sstring (glob_path));
		    }
		  globfree (&glob_str);
		}

	      // loop over glob expansions
	      for_const_contents (vector<sstring>, globbed_includes, arg)
		{
		  CTAG(2,SEXPR_MACROS SEXPR_INCLUDE) << "Including file " << *arg << '\n';
		  
		  SExpr_file included_file (arg->c_str());
		  parent_sexpr.child.insert (child_iter, included_file.sexpr.child.begin(), included_file.sexpr.child.end());
		}
	      erase_pos.push_back (child_iter);
	    }
	}
    }

  // erase the include statements
  for_contents (list<SExprIter>, erase_pos, erase_iter)
    parent_sexpr.child.erase (*erase_iter);
}


void SExpr_macros::visit (SExpr& parent_sexpr)
{
  // process substitutions at top level
  handle_replace (parent_sexpr);

  // look at children
  visit_and_reap (parent_sexpr);
}

void SExpr_macros::visit_and_reap (SExpr& parent_sexpr)
{
  list<SExprIter> erase;  // temporary list of positions to erase
  for_iterator (SExprIter, child_iter, parent_sexpr.child.begin(), parent_sexpr.child.end())
    visit_child (parent_sexpr, child_iter, erase);
  for_contents (list<SExprIter>, erase, erase_iter)
    parent_sexpr.child.erase (*erase_iter);
}

void SExpr_macros::handle_replace (SExpr& sexpr)
{
  if (sexpr.is_atom() && replace.find (sexpr.atom) != replace.end())
    sexpr.atom = replace[sexpr.atom];
}

void SExpr_macros::visit_child (SExpr& parent_sexpr, SExprIter& child_iter, list<SExprIter>& erase)
{
  SExpr& sexpr = *child_iter;
  // macros
  if (sexpr.is_list() && !sexpr.child.empty() && sexpr[0].is_atom())
    {
      // process substitutions one level down
      for_contents (list<SExpr>, sexpr.child, grandchild)
	handle_replace (*grandchild);
      // parse the operation
      sstring op (sexpr[0].atom);  // make a copy, not a reference, since we will be tampering with the expression
      CTAG(2,SEXPR_MACROS) << "Considering " << op << '\n';
      op = sexpr_macro_aliases.expand(op);
      if (op == SEXPR_WARN)       // warn
	{
	  // expand any macros in the warning text (disgusting hack)
	  SExpr_list_operations list_ops;
	  preorder_visit (sexpr);
	  list_ops.postorder_visit (sexpr);
	  // print the warning
	  CLOGERR << sexpr.atoms_to_strings(1) << '\n';
	  erase.push_back (child_iter);
	}
      else if (op == SEXPR_DEFINE)       // define
	{
	  CTAG(1,SEXPR_MACROS) << "Processing constant definition\n";
	  if (sexpr.child.size() == 3 && sexpr[1].is_atom() && sexpr[2].is_atom())
	    {
	      replace[sexpr[1].atom] = sexpr[2].atom;
	      erase.push_back (child_iter);
	    }
	  else
	    THROWEXPR ("Bad format for constant definition: " << sexpr);
	}
      else if (op == SEXPR_FOREACH)       // foreach
	{
	  CTAG(1,SEXPR_MACROS) << "Processing explicit iteration\n";
	  if (sexpr.child.size() >= 4 && sexpr[2].is_list())
	    {
	      // expand any macros in the arguments (hack hack hack)
	      SExpr_list_operations list_ops;
	      preorder_visit (sexpr[2]);
	      list_ops.postorder_visit (sexpr[2]);
	      // do the iteration
	      const vector<sstring> foreach_list = sexpr[2].atoms_to_strings();
	      expand_foreach (parent_sexpr, child_iter, 3, foreach_list, erase);
	    }
	  else
	    THROWEXPR ("Bad format for explicit iteration: " << sexpr);
	}
      else if (op == SEXPR_FOREACH_INT)       // foreach integer in range
	{
	  CTAG(1,SEXPR_MACROS) << "Processing iteration over integer range\n";
	  if (sexpr.child.size() >= 4 && sexpr[2].is_list() && sexpr[2].child.size() == 2 && sexpr[2][0].is_atom() && sexpr[2][1].is_atom())
	    {
	      // expand macros in the range expression (so hacky)
	      SExpr_list_operations list_ops;
	      for_contents (list<SExpr>, sexpr[2].child, great_grandchild)
		{
		  preorder_visit (*great_grandchild);
		  list_ops.postorder_visit (*great_grandchild);
		}
	      // do the iteration
	      vector<sstring> foreach_list;
	      for (int i = sexpr[2][0].atom.to_int(); i <= sexpr[2][1].atom.to_int(); ++i)
		{
		  sstring s;
		  s << i;
		  foreach_list.push_back (s);
		}
	      expand_foreach (parent_sexpr, child_iter, 3, foreach_list, erase);
	    }
	  else
	    THROWEXPR ("Bad format for integer range iteration: " << sexpr);
	}
      else if (foreach.find (op) != foreach.end())       // other iterations
	{	
	  CTAG(1,SEXPR_MACROS) << "Processing pre-defined iteration\n";
	  if (sexpr.child.size() >= 3)
	    expand_foreach (parent_sexpr, child_iter, 2, foreach[op], erase);
	  else
	    THROWEXPR ("Bad format for preset iteration: " << sexpr);
	}
    }
}

void SExpr_macros::expand_foreach (SExpr& parent_sexpr, SExprIter& parent_pos, unsigned int element_offset, const vector<sstring>& foreach_list, list<SExprIter>& erase)
{
  SExpr& foreach_sexpr = *parent_pos;
  if (CTAGGING(-2,SEXPR_MACROS))
    CL << "About to expand iteration: " << foreach_sexpr << '\n';

  if (!foreach_sexpr[1].is_atom())
    THROWEXPR ("Term following '" << foreach_sexpr[0].atom << "' should be a variable name, but is a list: " << foreach_sexpr);
  const sstring& foreach_var = foreach_sexpr[1].atom;
  if (replace.find (foreach_var) != replace.end())
    THROWEXPR ("Iterator variable '" << foreach_var << "' is already in use: " << foreach_sexpr);

  list<SExpr>::const_iterator
    element_begin = foreach_sexpr.child.begin(),
    element_end = foreach_sexpr.child.end();
  while (element_offset--)
    ++element_begin;
  for_const_contents (vector<sstring>, foreach_list, foreach_val)
    {
      CTAG(0,SEXPR_MACROS SEXPR_EXPAND) << "Expanding " << foreach_var << " to " << *foreach_val << '\n';
      replace[foreach_var] = *foreach_val;
      for_iterator (list<SExpr>::const_iterator, element_iter, element_begin, element_end)
	{
	  SExprIter new_child = parent_sexpr.child.insert (parent_pos, *element_iter);
	  const unsigned int old_erase_sz = erase.size();  // hacky way of finding out if erase has changed
	  visit_child (parent_sexpr, new_child, erase);
	  // if visit_child added new_child to erase(), then we shouldn't try to visit it.
	  // yeah, this is so hacky and brittle it may not even survive you looking at it. write a proper lisp interpreter.
	  if (old_erase_sz == erase.size())
	    preorder_visit (*new_child);
	}
    }
  replace.erase (foreach_var);
  erase.push_back (parent_pos);
}

void SExpr_list_operations::visit (SExpr& parent_sexpr)
{
  for_iterator (list<SExpr>::iterator, child_iter, parent_sexpr.child.begin(), parent_sexpr.child.end())
    {
      SExpr& sexpr = *child_iter;
      if (sexpr.is_list() && !sexpr.child.empty() && sexpr[0].is_atom())
	{
	  sstring op (sexpr[0].atom);  // make a copy, not a reference, since we will be tampering with the expression
	  CTAG(2,SEXPR_MACROS) << "Considering " << op << '\n';
	  op = sexpr_macro_aliases.expand(op);
	  // list operations
	  if (op == SEXPR_CONCATENATE || op == SEXPR_SUM || op == SEXPR_MULTIPLY || op == SEXPR_MODULUS
	      || op == SEXPR_AND || op == SEXPR_OR)
	    {
	      CTAG(1,SEXPR_MACROS) << "Processing list operation\n";
	      vector<sstring> atoms = sexpr.atoms_to_strings (1);
	      // wipe the list expression
	      sexpr.child.clear();
	      sexpr.atom.clear();
	      // perform the list operation
	      if (op == SEXPR_CONCATENATE)
		sexpr.atom = sstring::join (atoms, "");
	      else if (op == SEXPR_SUM)
		{
		  double tot = 0.;
		  for_const_contents (vector<sstring>, atoms, a)
		    tot += a->to_double();
		  sexpr.atom << tot;
		}
	      else if (op == SEXPR_MULTIPLY)
		{
		  double tot = 1.;
		  for_const_contents (vector<sstring>, atoms, a)
		    tot *= a->to_double();
		  sexpr.atom << tot;
		}
	      else if (op == SEXPR_AND)
		{
		  bool result = true;
		  for_const_contents (vector<sstring>, atoms, a)
		    result = result && a->to_int();
		  sexpr.atom = result ? SEXPR_TRUE : SEXPR_FALSE;
		}
	      else if (op == SEXPR_OR)
		{
		  bool result = false;
		  for_const_contents (vector<sstring>, atoms, a)
		    result = result || a->to_int();
		  sexpr.atom = result ? SEXPR_TRUE : SEXPR_FALSE;
		}
	      else if (op == SEXPR_MODULUS)
		{
		  if (atoms.size() != 2)
		    THROWEXPR ("Modulus needs exactly two arguments");
		  sexpr.atom << (atoms[0].to_int() % atoms[1].to_int());
		}
	      else
		THROWEXPR ("Should never get here");
	    }
	  // functions with one atomic argument
	  else if (op == SEXPR_INT || op == SEXPR_CHR || op == SEXPR_ORD || op == SEXPR_NOT)
	    {
	      CTAG(1,SEXPR_MACROS) << "Processing " << op << '\n';
	      if (sexpr.child.size() != 2 || !sexpr[1].is_atom())
		THROWEXPR ("Expected one integer argument for " << op << ": " << sexpr);
	      const sstring arg = sexpr[1].atom;
	      sexpr.child.clear();
	      sexpr.atom.clear();
	      if (op == SEXPR_INT)
		sexpr.atom << arg.to_int();
	      else if (op == SEXPR_CHR)
		sexpr.atom << (char) arg.to_int();
	      else if (op == SEXPR_ORD)
		sexpr.atom << (int) arg[0];
	      else if (op == SEXPR_NOT)
		sexpr.atom << (arg.to_int() ? SEXPR_FALSE : SEXPR_TRUE);
	      else
		THROWEXPR ("Should never get here");
	    }
	  // functions with two atomic arguments
	  else if (op == SEXPR_DIVIDE || op == SEXPR_SUBTRACT || op == SEXPR_GREATER || op == SEXPR_LESS || op == SEXPR_GEQ || op == SEXPR_LEQ)
	    {
	      CTAG(1,SEXPR_MACROS) << "Processing " << op << '\n';
	      if (sexpr.child.size() != 3 || !sexpr[1].is_atom() || !sexpr[2].is_atom())
		THROWEXPR ("Expected two arguments for " << op << " function: " << sexpr);
	      const sstring arg1 = sexpr[1].atom, arg2 = sexpr[2].atom;
	      sexpr.child.clear();
	      sexpr.atom.clear();
	      if (op == SEXPR_DIVIDE)
		sexpr.atom << (double) (arg1.to_double() / arg2.to_double());
	      else if (op == SEXPR_SUBTRACT)
		sexpr.atom << (double) (arg1.to_double() - arg2.to_double());
	      else if (op == SEXPR_GREATER)
		sexpr.atom << (arg1.to_double() > arg2.to_double() ? SEXPR_TRUE : SEXPR_FALSE);
	      else if (op == SEXPR_LESS)
		sexpr.atom << (arg1.to_double() < arg2.to_double() ? SEXPR_TRUE : SEXPR_FALSE);
	      else if (op == SEXPR_GEQ)
		sexpr.atom << (arg1.to_double() >= arg2.to_double() ? SEXPR_TRUE : SEXPR_FALSE);
	      else if (op == SEXPR_LEQ)
		sexpr.atom << (arg1.to_double() <= arg2.to_double() ? SEXPR_TRUE : SEXPR_FALSE);
	      else
		THROWEXPR ("Should never get here");
	    }
	  // functions with two S-expression arguments returning true or false
	  else if (op == SEXPR_EQUALS || op == SEXPR_NOT_EQUALS)
	    {
	      CTAG(1,SEXPR_MACROS) << "Processing " << op << '\n';
	      if (sexpr.child.size() != 3)
		THROWEXPR ("Expected three terms in " << op << " function: " << sexpr);
	      bool result;
	      if (op == SEXPR_EQUALS)
		result = sexpr[1] == sexpr[2];
	      else if (op == SEXPR_NOT_EQUALS)
		result = sexpr[1] != sexpr[2];
	      else
		THROWEXPR ("Should never get here");
	      sexpr.child.clear();
	      sexpr.atom.clear();
	      sexpr.atom << (result ? SEXPR_TRUE : SEXPR_FALSE);
	    }
	  // conditional function: three arguments, incl. two atoms + one S-expression
	  else if (op == SEXPR_CONDITIONAL)
	    {
	      CTAG(1,SEXPR_MACROS) << "Processing conditional\n";
	      if (sexpr.child.size() != 4 || !sexpr[1].is_atom())
		THROWEXPR ("Bad format for conditional: " << sexpr);
	      SExpr result_sexpr;
	      result_sexpr.swap (sexpr[1].get_atom().to_int() ? sexpr[2] : sexpr[3]);
	      sexpr.swap (result_sexpr);
	    }
	}
    }
}

#ifdef GUILE_INCLUDED
// dummy function to pass to scm_with_guile
static void*
register_functions (void* data)
{
  return NULL;
}
#endif /* GUILE_INCLUDED */

// SExpr_Scheme_evaluator
SExpr_Scheme_evaluator::SExpr_Scheme_evaluator()
{
#ifdef GUILE_INCLUDED
  scm_with_guile (&register_functions, NULL);
  write_proc = scm_variable_ref (scm_c_lookup("write"));
#endif /* GUILE_INCLUDED */
}

void SExpr_Scheme_evaluator::expand_Scheme_expressions (SExpr& sexpr) const
{
  typedef SExpr::SExprIter SExprIter;
  list<SExprIter> erase_pos;
  for_contents (list<SExpr>, sexpr.child, c)
    {
      bool is_scheme = false;
      if (c->is_list() && !c->child.empty() && (*c)[0].is_atom())
	{
	  const SExpr_atom& op ((*c)[0].get_atom());
	  if (op == SEXPR_EVAL)
	    {
	      // evaluate &scheme arguments as Scheme expression, splice results into sexpr.child before c
	      SExpr result = evaluate(*c);
	      if (result.is_atom())
		c->child.insert (c, result);
	      else
		c->child.insert (c, result.child.begin(), result.child.end());
	      is_scheme = true;
	    }
	  else if (op == SEXPR_EXEC)
	    {
	      // evaluate &exec arguments as Scheme expressions, discard results
	      evaluate(*c);
	      is_scheme = true;
	    }
	}
      if (is_scheme)
	erase_pos.push_back (c);  // mark block for deletion
      else
	expand_Scheme_expressions (*c);  // recursively descend
    }

  // erase the SExpr's marked for deletion
  for_contents (list<SExprIter>, erase_pos, erase_iter)
    sexpr.child.erase (*erase_iter);
}

SExpr SExpr_Scheme_evaluator::evaluate (SExpr& sexpr) const
{
  if (sexpr.is_atom())
    THROWEXPR("SExpr_Scheme_evaluator::evaluate called on an atom");
  if (sexpr.is_empty_list())
    THROWEXPR("SExpr_Scheme_evaluator::evaluate called on an empty list");
#ifdef GUILE_INCLUDED
  sstring result_string;
  // iterate over all but the first child
  SExpr::SExprIter iter = sexpr.child.begin();
  for (++iter; iter != sexpr.child.end(); ++iter)
    {
      sstring sexpr_string = iter->to_parenthesized_string();
      // evaluate; convert to C string
      CTAG(3,GUILE) << "Evaluating Scheme expression: " << sexpr_string << '\n';
      SCM result_scm = scm_c_eval_string (sexpr_string.c_str());
      SCM result_scm_string = scm_object_to_string (result_scm, write_proc);
      const char* result_c_string = scm_to_locale_string (result_scm_string);
      CTAG(3,GUILE) << "Result of Scheme evaluation: " << result_c_string << '\n';
      // catch SCM_UNSPECIFIED
      if (result_scm == SCM_UNSPECIFIED)
	CTAG(3,GUILE) << "Discarding unspecified return value from Scheme evaluation\n";
      else
	result_string << result_c_string << ' ';
    }
  // convert to SExpr, and return
  const SExpr result_sexpr (result_string.begin(), result_string.end());
  return result_sexpr;
#else /* GUILE_INCLUDED */
  THROWEXPR("This program was not compiled with the GNU Guile library, so " << SEXPR_EVAL << " and " << SEXPR_EXEC << " are not available.\n"
	    << "Recompile with Guile enabled to access these Scheme functions.");
  return SExpr();
#endif /* GUILE_INCLUDED */
}

