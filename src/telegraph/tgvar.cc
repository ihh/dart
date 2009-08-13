#include <fstream>
#include <sstream>

#include "telegraph/tgvar.h"
#include "util/Regexp.h"

Telegraph_PScores_adaptor::Telegraph_PScores_adaptor (PScores& pscores)
  : pscores (pscores)
{
  for (int g = 0; g < pscores.pgroups(); ++g)
    io_groups.insert (g);
}

Telegraph_PScores_adaptor::Telegraph_PScores_adaptor (Dirichlet_prior& prior)
  : pscores (*prior.pscores)
{
  for_const_contents (Dirichlet_prior::Assignment_set, prior.assignment, multi_ass)
    for_const_contents (Multigroup, multi_ass->first, g)
    io_groups.insert (g->group_idx);
}

sstring Telegraph_PScores_adaptor::telegraph_varname (int g, int v) const
{
  sstring varname;
  varname << pscores.group_name[g] << '[';
  if ((int) pscores.group_suffix[g].size() > v)
    varname << pscores.group_suffix[g][v];
  else
    varname << v;
  varname << ']';
  return varname;
}

// regexp to recognise "GROUP[VAR] => VALUE"
Regexp assign_re = "^[ \t\n]*([A-Za-z_][A-Za-z_0-9]*)[ \t\n]*\\[[ \t\n]*([A-Za-z_0-9]+)[ \t\n]*\\][ \t\n]*=>[ \t\n]*(.+)$";

// regexp to recognise "(int).(int)"
Regexp bits_re = "^([\\+\\-]?[0-9]*)[ \t\n]*\\.[ \t\n]*([0-9]+)[ \t\n]*$";

void Telegraph_PScores_adaptor::read (const char* filename, bool ignore_unknown_params, bool ignore_undefined_params)
{
  CTAG(5,TELEGRAPH) << "Reading parameters from Telegraph file '" << filename << "'\n";
  ifstream in (filename);
  if (!in) THROWEXPR ("Telegraph parameter file '" << filename << "' not found\n");
  read (in, ignore_unknown_params, ignore_undefined_params);
}

void Telegraph_PScores_adaptor::read_from_string (const char* param_str, bool ignore_unknown_params, bool ignore_undefined_params)
{
  const sstring s (param_str);
  istringstream in (s);
  read (in, ignore_unknown_params, ignore_undefined_params);
}

void Telegraph_PScores_adaptor::read (istream& in, bool ignore_unknown_params, bool ignore_undefined_params)
{
  // read in entire file as a string
  sstring instr, s;
  while (in && !in.eof())
    {
      s.getline (in);
      s.chomp();
      instr << s;
    }
  // split into individual statements at semicolons
  const vector<sstring> lines = instr.split (";");
  sstring groupname_string, varname_string, gvname_string, value_string;
  sstring div_string, mod_string;
  // figure out every legal VARNAME string from the PScores object (wasteful but easy)
  typedef map<sstring,PVar> string2pvar_map;
  string2pvar_map str2pvar;
  for_const_contents (set<int>, io_groups, g)
    for (int v = 0; v < pscores.pgroup_size(*g); ++v)
      str2pvar [telegraph_varname(*g,v)] = PVar(*g,v);
  set<PVar> seen_pvar;
  // loop through all statements
  for_const_contents (vector<sstring>, lines, line)
    {
      CTAG(2,TELEGRAPH_PARSER) << "Parsing line '" << *line << "'\n";
      if (assign_re.Match (line->c_str()))
	{
	  // extract GROUPNAME, VARNAME and VALUE strings
	  groupname_string = assign_re[1];
	  varname_string = assign_re[2];
	  value_string = assign_re[3];
	  // make string GVNAME = "GROUP[VAR]"
	  gvname_string.clear();
	  gvname_string << groupname_string << '[' << varname_string << ']';
	  // map the GVNAME to one of DART's PVar's
	  if (str2pvar.find (gvname_string) == str2pvar.end())
	    {
	      if (ignore_unknown_params)
		continue;
	      else
		THROWEXPR ("Unknown parameter name in Telegraph parameter file: " << gvname_string);
	    }
	  const PVar pv = str2pvar[gvname_string];
	  // check we haven't already seen this PVar
	  if (seen_pvar.find (pv) != seen_pvar.end())
	    THROWEXPR ("Duplicate parameter assignment in Telegraph parameter file: " << gvname_string);
	  seen_pvar.insert (pv);
	  // parse bitscore format
	  if (bits_re.Match (value_string.c_str()))
	    {
	      div_string = bits_re[1];
	      mod_string = bits_re[2];
	      const int sign = div_string[0] == '-' ? +1 : -1;
	      pscores[pv] = sign * (abs(div_string.to_int()) * DartScore2BitsRatio + mod_string.to_int());
	    }
	}
      else
	THROWEXPR ("Bad statement format in Telegraph parameter file: " << *line);
    }
  // check if any parameters remain unassigned; if so, squeal
  if (!ignore_undefined_params)
    if (seen_pvar.size() != str2pvar.size())
      {
	sstring unseen;
	for_const_contents (string2pvar_map, str2pvar, sp)
	  if (seen_pvar.find (sp->second) == seen_pvar.end())
	    unseen << ' ' << sp->first;
	CLOGERR << "Warning: Telegraph parameter file left the following parameters unassigned:" << unseen << "\n";
      }
}

void Telegraph_PScores_adaptor::write (const char* filename) const
{
  ofstream out (filename);
  write (out);
}

void Telegraph_PScores_adaptor::write (ostream& out) const
{
  // figure out the precision
  sstring tmp;
  tmp << (int) (DartScore2BitsRatio-1);
  const int prec = tmp.size();
  // output the vars
  for_const_contents (set<int>, io_groups, g)
    {
      const int group_size = pscores.pgroup_size(*g);
      for (int v = 0; v < group_size; ++v)
	{
	  const Score sc = pscores.group[*g][v];
	  const int div = abs(sc) / DartScore2BitsRatio;
	  const int mod = abs(sc) % DartScore2BitsRatio;
	  out << telegraph_varname(*g,v) << " => ";
	  if (sc > 0) out << '-';
	  out << div << '.';
	  out.fill('0');
	  out.width(prec);
	  out << mod << ';';
	  out << (v < group_size-1 ? '\t' : '\n');
	}
    }
}
