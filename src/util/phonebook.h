#ifndef PHONEBOOK_INCLUDED
#define PHONEBOOK_INCLUDED

#include <map>
#include "util/macros.h"

// Phonebook class (map from sstring to int)

class Phonebook : public map<sstring,int>
{
public:
  template <class InputIterator>
  void make_index (InputIterator first, InputIterator last, int first_n = 0)
    { for_iterator (InputIterator, i, first, last) insert (value_type ((*i), first_n++)); }

  Phonebook() : map<sstring,int>() { }

  template <class InputIterator>
  Phonebook (InputIterator first, InputIterator last, int first_n = 0) : map<sstring,int>() { make_index (first, last, first_n); }
  
  // contains() method
  
  bool contains (const sstring& k) const
    { const_iterator i = find (k); return i != end(); }

  // lookup() method
  
  const int lookup (const sstring& k, const int default_value = -1) const
    { const_iterator i = find (k); return i == end() ? default_value : (*i).second; }
};

#endif
