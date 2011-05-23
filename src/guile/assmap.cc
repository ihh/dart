#include "guile/assmap.h"
#include "util/guile-defs.h"
#include "util/logfile.h"

Ass_map::Ass_map (SCM scm_list)
{
  /* Check that first argument is a list */ 
  if (SCM_NFALSEP(scm_list_p(scm_list)))
    {
      /* Iterate through the list */ 
      while (SCM_FALSEP (scm_null_p (scm_list)))
	{
	  /* Get the head of the list */ 
	  SCM head = SCM_CAR(scm_list);
	  // Be forgiving of dart's SExpr syntax, which doesn't distinguish between pairs and lists:
	  // allow either a pair, or a two-element list
	  SCM tag, value;
	  if (SCM_NFALSEP(scm_list_p(head)))
	    {
	      tag = SCM_CAR(head);
	      value = SCM_CADR(head);
	    }
	  else
	    if (SCM_NFALSEP(scm_pair_p(head)))
	      {
		tag = SCM_CAR(head);
		value = SCM_CDR(head);
	      }
	    else
	      CLOGERR << "Discarding " << scm_to_string(head) << "\n";

	  // get the string
	  if (!scm_is_string(tag))
	    CTAG(3,GUILE) << "In association-list wrapper: Scheme expression is not a string\n";

	  const sstring tag_str = scm_to_string (tag);

	  /* Add it to the map */
	  (*this) [tag_str] = value;

	  /* Discard the head of the list */
	  scm_list = SCM_CDR(scm_list);
	}
    }
  else   // first arg is not a list
    THROWEXPR ("Attempt to construct association-list wrapper from non-list SCM object");
}

Ass_map::Ass_map (SExpr& sexpr, int offset)
{
  if (sexpr.is_list())
    {
      int n = 0;
      for_contents (list<SExpr>, sexpr.child, child)
	{
	  if (++n >= offset)
	    {
	      if (child->has_tag() && child->has_value())
		(*this) [child->tag()] = sexpr_to_scm (&child->value());
	      else
		CLOGERR << "Discarding " << child->to_parenthesized_string() << "\n";
	    }
	}
    }
  else
    THROWEXPR ("Attempt to construct association-list wrapper from non-list SExpr object");
}

SCM Ass_map::value_or_false (const char* tag) const
{
  const sstring tag_str (tag);
  const const_iterator iter = find (tag_str);
  if (iter == end())
    return SCM_BOOL_F;
  return iter->second;
}

SCM Ass_map::value_or_die (const char* tag) const
{
  const sstring tag_str (tag);
  const const_iterator iter = find (tag_str);
  if (iter == end())
    THROWEXPR ("In value_or_die: Couldn't find " << tag);
  return iter->second;
}
