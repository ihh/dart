#ifndef MAP_KEYS_INCLUDED
#define MAP_KEYS_INCLUDED

#include <map>
#include <vector>
#include <algorithm>

template <class Key, class T, class Compare, class Alloc>
vector<Key> map_keys (const map<Key,T,Compare,Alloc>& m)
{
  typedef map<Key,T,Compare,Alloc> MapType;
  vector<Key> keys;
  keys.reserve (m.size());
  template_for_const_contents (MapType, m, keyval)
    keys.push_back (keyval->first);
  return keys;
}

template <class Key, class T, class Compare, class Alloc>
vector<T> map_values (const map<Key,T,Compare,Alloc>& m)
{
  typedef map<Key,T,Compare,Alloc> MapType;
  vector<T> values;
  values.reserve (m.size());
  template_for_const_contents (MapType, m, keyval)
    values.push_back (keyval->second);
  return values;
}

#endif
