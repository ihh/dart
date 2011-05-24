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

	  // get the tag string
	  sstring tag_str;
	  if (scm_is_string(tag))
	    tag_str = scm_to_string_unquoted (tag);
	  else if (scm_is_symbol(tag))
	    tag_str = scm_to_string_unquoted (scm_symbol_to_string (tag));
	  else
	    CTAG(3,GUILE) << "In association-list wrapper: Scheme expression " << scm_to_string(tag) << " is not a string or a symbol\n";

	  /* Add it to the map */
	  if (tag_str.size())
	    {
	      CTAG(1,GUILE) << "Parsed parameter " << tag_str << '\n';
	      (*this) [tag_str] = value;
	    }

	  /* Discard the head of the list */
	  scm_list = SCM_CDR(scm_list);
	}
    }
  else
    {
      if (SCM_NFALSEP(scm_hash_table_p(scm_list)))
	{
	  THROWEXPR ("Attempt to construct association-list wrapper from a hash-table SCM object; this is currently not supported");
	}
      else   // first arg is not a recognized type
	THROWEXPR ("Attempt to construct association-list wrapper from unrecognized non-list SCM object");
    }
}

Ass_map::Ass_map (SExpr& sexpr, int offset)
{
  if (sexpr.is_list())
    {
      int n = 0;
      for_contents (list<SExpr>, sexpr.child, child)
	{
	  if (++n > offset)
	    {
	      // be forgiving: accept (tag value) or (tag . value)
	      if (child->has_tag() && child->has_value())  // looks like a (tag value) list?
		(*this) [child->tag()] = sexpr_to_scm (&child->value());
	      else if (child->has_tag() && child->is_scheme_pair())  // looks like a Scheme (tag . value) pair?
		(*this) [child->tag()] = sexpr_to_scm (&((*child)[2]));
	      else
		CLOGERR << "Discarding " << child->to_parenthesized_string() << "\n";
	    }
	}
    }
  else
    THROWEXPR ("Attempt to construct association-list wrapper from non-list SExpr object");
}

Ass_map::const_iterator Ass_map::find (const char* tag) const
{
  const sstring tag_str (tag);
  const const_iterator iter = map<sstring,SCM>::find (tag_str);
  return iter;
}

SCM Ass_map::scm_value_or_false (const char* tag) const
{
  const const_iterator iter = find (tag);
  if (iter == end())
    return SCM_BOOL_F;
  return iter->second;
}

SCM Ass_map::scm_value_or_die (const char* tag) const
{
  const const_iterator iter = find (tag);
  if (iter == end())
    THROWEXPR ("In scm_value_or_die: Couldn't find " << tag);
  return iter->second;
}

SExpr Ass_map::sexpr_value_or_empty_list (const char* tag) const
{
  const const_iterator iter = find (tag);
  if (iter == end())
    return SExpr();
  return scm_to_sexpr (iter->second);
}

SExpr Ass_map::sexpr_value_or_die (const char* tag) const
{
  return scm_to_sexpr (scm_value_or_die (tag));
}

SExpr* Ass_map::new_sexpr_value_or_null (const char* tag) const
{
  const const_iterator iter = find (tag);
  if (iter == end())
    return NULL;
  return scm_to_new_sexpr (iter->second);
}
