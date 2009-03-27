#ifndef DART_MAP_INCLUDED
#define DART_MAP_INCLUDED

#include <map>
#include "util/map_keys.h"

// map class with contains() and lookup()

#ifndef __STL_LIMITED_DEFAULT_TEMPLATES
template <class Key, class T, class Compare = less<Key>, class Alloc = alloc>
#else
template <class Key, class T, class Compare, class Alloc = alloc>
#endif
class dartmap : public map <Key, T, Compare, Alloc>
{
public:
  
  // wrappers for map constructors
  
  ymap() : map() {}
  explicit ymap(const Compare& comp) : map (comp) {}

#ifdef __STL_MEMBER_TEMPLATES
  template <class InputIterator>
  ymap(InputIterator first, InputIterator last)
    : map (first, last) { }

  template <class InputIterator>
  ymap(InputIterator first, InputIterator last, const Compare& comp)
    : map (first, last, comp) { }
#else
  ymap(const value_type* first, const value_type* last)
    : map (first, last) { }
  ymap(const value_type* first, const value_type* last, const Compare& comp)
    : map (first, last, comp) { }

  ymap(const_iterator first, const_iterator last)
    : map (first, last) { }
  ymap(const_iterator first, const_iterator last, const Compare& comp)
    : map (first, last, comp) { }
#endif /* __STL_MEMBER_TEMPLATES */

  // contains() method
  
  bool contains (const key_type& k) const
    { const_iterator i = find (k); return i != end(); }

  // lookup() method
  
  const data_type& lookup (const key_type& k, const data_type& default_value) const
    { const_iterator i = find (k); return i == end() ? default_value : (*i).second; }

  // keys() & values() methods

  vector<key_type>   keys() const { return map_keys (*this); }
  vector<value_type> values() const { return map_values (*this); }
};

#endif
