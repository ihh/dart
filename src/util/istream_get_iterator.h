// istream_iterator that uses get(char&) instead of operator>>(char&)

#ifndef ISTREAM_GET_ITERATOR
#define ISTREAM_GET_ITERATOR

class istream_get_iterator
{
  friend inline bool operator==(const istream_get_iterator&__x,
				const istream_get_iterator&__y) {
    return (__x._M_stream == __y._M_stream &&
	    __x._M_end_marker == __y._M_end_marker) ||
      __x._M_end_marker == false && __y._M_end_marker == false;
  }
  
  friend inline bool operator!=(const istream_get_iterator&__x,
				const istream_get_iterator&__y)
  { return ! (__x == __y); }

 protected:
  istream* _M_stream;
  char _M_value;
  bool _M_end_marker;

  void _M_read() {
    _M_end_marker = (*_M_stream) ? true : false;
    if (_M_end_marker) _M_stream->get(_M_value);
    _M_end_marker = (*_M_stream) ? true : false;
  }

 public:
  typedef input_iterator_tag             iterator_category;
  typedef char                           value_type;
  typedef ptrdiff_t                      difference_type;
  typedef const char*                    pointer;
  typedef const char&                    reference;

  istream_get_iterator() : _M_stream(&cin), _M_end_marker(false) {}

  istream_get_iterator(istream& __s) : _M_stream(&__s) { _M_read(); }
  reference operator*() const { return _M_value; }
#ifndef __SGI_STL_NO_ARROW_OPERATOR
  pointer operator->() const { return &(operator*()); }
#endif /* __SGI_STL_NO_ARROW_OPERATOR */
  istream_get_iterator& operator++() { 
    _M_read(); 
    return *this;
  }
  istream_get_iterator operator++(int)  {
    istream_get_iterator __tmp = *this;
    _M_read();
    return __tmp;
  }

};




#endif

