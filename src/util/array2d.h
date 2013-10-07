#ifndef ARRAY2D_INCLUDED
#define ARRAY2D_INCLUDED

#include <vector>
#include <functional>
#include <algorithm>
#include "util/dexception.h"
#include "util/logfile.h"

using namespace std;

// xy_container's must implement the following
//   operator==
//   swap (xy_container_type)
//   default, copy constructors
//  ...plus the following xy constructor...
//   xy_container_type (int xsize, int ysize, const_reference default_val)
//  ...and the following xy accessors...
//   inline const_reference get_xy_or_default (int x, int y, int xsize, int ysize) const    [nonperturbative]
//   inline reference get_xy (int x, int y, int xsize, int ysize)  [creates (x,y) entry if it does not exist]

template<class T>
struct array2d_dense_vector : vector<T>
{
  inline array2d_dense_vector() : vector<T>() { }
  inline array2d_dense_vector (const array2d_dense_vector<T>& a) : vector<T> (a) { }
  inline array2d_dense_vector (int xsize, int ysize, const T& default_val) : vector<T> (xsize*ysize, default_val) { }
  inline const T& get_xy_or_default (int x, int y, int xsize, int ysize) const { return (*this) [x + y * xsize]; }
  inline T& get_xy (int x, int y, int xsize, int ysize) { return (*this) [x + y * xsize]; }
};

template<class T>
struct array2d_sparse_vector : vector< map<int,T> >
{
  // typedefs
  typedef typename map<int,T>::iterator map_iterator;
  typedef typename map<int,T>::const_iterator const_map_iterator;

  // data: default value
  T default_val;

  // methods
  inline array2d_sparse_vector() : vector< map<int,T> >() { }
  inline array2d_sparse_vector (const array2d_sparse_vector<T>& a) : vector< map<int,T> > (a) { }
  inline array2d_sparse_vector (int xsize, int ysize, const T& default_val)
    : vector< map<int,T> > (ysize),
      default_val (default_val)
  { }

  inline const T& get_xy_or_default (int x, int y, int xsize, int ysize) const {
    const map<int,T>& col = (*this)[y];
    const_map_iterator i = col.find (x);
    if (i == col.end())
      return default_val;
    return (*i).second;
  }

  inline T& get_xy (int x, int y, int xsize, int ysize) {
    map<int,T>& col = (*this)[y];
    map_iterator i = col.find (x);
    if (i == col.end())
      return col[x] = default_val;
    return (*i).second;
  }
};

// array2d class

/* Problems: confusingly, operator<< and operator>> output the transpose of the array, i.e. rows become columns
   This confusion arises from the (x,y) indexing convention of arrays and the (row,col) indexing of matrices.
   As a workaround, use write_rowcol() and read_rowcol() methods.
*/

template <class T, class xy_container_type = array2d_dense_vector<T> >
class array2d
{
public:
  typedef T          value_type;
  typedef T&         reference;
  typedef const T&   const_reference;
  typedef T*         pointer;
  
  typedef ptrdiff_t  difference_type;
  typedef size_t     size_type;

  typedef array2d <value_type, xy_container_type> array2d_type;

 private:
  int _xsize, _ysize;
  xy_container_type _data;
  
 public:
  
  vector<int> dim() const { vector<int> d(2); d[0] = _xsize; d[1] = _ysize; return d; }
  int xsize() const { return _xsize; }
  int ysize() const { return _ysize; }
  int columns() const { return _xsize; }
  int rows() const { return _ysize; }
  int rank() const { return 2; }

  bool operator== (const array2d_type& a) const
  { return _xsize == a._xsize && _ysize == a._ysize && _data == a._data; }

  // NB row/column addressing is opposite to matrices:
  // matrices use (row,column), we use (x,y)=(column,row)
  //

  reference operator()(int x, int y)
  {
#ifdef DART_DEBUG
    if (x < 0 || x >= _xsize || y < 0 || y >= _ysize)
      DART_DEBUG_ERROR ("array2d overflow");
#endif
    return _data.get_xy (x, y, _xsize, _ysize);
  }
  const_reference operator()(int x, int y) const
  {
#ifdef DART_DEBUG
    if (x < 0 || x >= _xsize || y < 0 || y >= _ysize)
      DART_DEBUG_ERROR ("array2d overflow");
#endif
    return _data.get_xy_or_default (x, y, _xsize, _ysize);
  }

  reference operator[](const vector<int>& k) { return (*this) (k[0], k[1]); }
  const_reference operator[](const vector<int>& k) const { return (*this) (k[0], k[1]); }

  reference operator[](const pair<int, int>& k) { return (*this) (k.first, k.second); }
  const_reference operator[](const pair<int, int>& k) const { return (*this) (k.first, k.second); }

  reference entry (int x, int y) { return (*this) (x, y); }
  const_reference entry (int x, int y) const { return (*this) (x, y); }

  void fill (const value_type& fill_value)
  {
    for (int x = 0; x < xsize(); ++x)
      for (int y = 0; y < ysize(); ++y)
	(*this)(x,y) = fill_value;
  }

  void resize (int new_xsize, int new_ysize, const value_type& default_value)
    {
      array2d_type new_array (new_xsize, new_ysize, default_value);
      for (int x = 0; x < min (_xsize, new_xsize); ++x)
	for (int y = 0; y < min (_ysize, new_ysize); ++y)
	  new_array (x, y) = entry (x, y);
      swap (new_array);
    }

  void swap (array2d_type& a)
    {
      std::swap (_xsize, a._xsize);
      std::swap (_ysize, a._ysize);
      _data.swap (a._data);
    }

  // iterators
  struct xy_iterator : public iterator <forward_iterator_tag, T>
  {
    // data
    int x, y;
    array2d_type* a;
    // methods
    xy_iterator() : x(0), y(0), a(0) { }
    xy_iterator (int x, int y, array2d_type* a) : x(x), y(y), a(a) { }
    T& operator*() { return (*a)(x,y); }
    const T& operator*() const { return (*a)(x,y); }
    bool operator== (const xy_iterator& i) const { return y==i.y && x==i.x && a==i.a; }
    bool operator!= (const xy_iterator& i) const { return !(*this==i); }
    xy_iterator& operator++() { if (++x >= a->xsize()) { x = 0; ++y; } return *this; }
  };

  typedef xy_iterator iterator;
  typedef xy_iterator const_iterator;

  iterator begin() { return xy_iterator (0, 0, this); }
  iterator end() { return xy_iterator (0, ysize(), this); }
  const_iterator begin() const { return xy_iterator (0, 0, (array2d_type*) this); }
  const_iterator end() const { return xy_iterator (0, ysize(), (array2d_type*) this); }

  // row & column iterators
  class Row_iterator
  {
    // iterator typedefs (inheriting from iterator seems to crash the compiler)
  public:
    typedef random_access_iterator_tag iterator_category;
    typedef T value_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;

    // the real stuff
  private:
    array2d_type* _array;
    int _row;
    int _col_offset;
    bool _same_row (const Row_iterator& i) const { return _array==i._array && _row==i._row; }
  public:
    Row_iterator (array2d_type& a, int row, int col_offset = 0) : _array(&a), _row(row), _col_offset(col_offset) { }
    reference operator[] (int col) { return (*_array) (_row, col + _col_offset); }
    const_reference operator[] (int col) const { return (*_array) (_row, col + _col_offset); }
    reference operator*() { return (*this)[0]; }
    const_reference operator*() const { return (*this)[0]; }
    bool operator== (const Row_iterator& i) const { return _same_row(i) && _col_offset==i._col_offset; }
    Row_iterator& operator++() { ++_col_offset; return *this; }
    Row_iterator operator++ (int) { Row_iterator tmp = *this; ++_col_offset; return tmp; }
    Row_iterator& operator--() { --_col_offset; return *this; }
    Row_iterator operator-- (int) { Row_iterator tmp = *this; --_col_offset; return tmp; }
    Row_iterator& operator+= (ptrdiff_t d) { _col_offset += d; return *this; }
    Row_iterator& operator-= (ptrdiff_t d) { _col_offset -= d; return *this; }
    friend Row_iterator operator+ (const Row_iterator& l, ptrdiff_t d) { Row_iterator res = l; res += d; return res; }
    friend Row_iterator operator- (const Row_iterator& l, ptrdiff_t d) { Row_iterator res = l; res -= d; return res; }
    ptrdiff_t operator- (const Row_iterator& i) const
    {
      if (!_same_row(i)) THROW Standard_exception ("Can't subtract array2d<>::Row_iterator's of different types");
      return _col_offset - i._col_offset;
    }
  };

  class Column_iterator
  {
    // iterator typedefs (inheriting from iterator seems to crash the compiler)
  public:
    typedef random_access_iterator_tag iterator_category;
    typedef T value_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;

    // the real stuff
  private:
    array2d_type* _array;
    int _col;
    int _row_offset;
    bool _same_col (const Column_iterator& i) const { return _array==i._array && _col==i._col; }
  public:
    Column_iterator (array2d_type& a, int col, int row_offset = 0) : _array(&a), _col(col), _row_offset(row_offset) { }
    reference operator[] (int row) { return (*_array) (row + _row_offset, _col); }
    const_reference operator[] (int row) const { return (*_array) (row + _row_offset, _col); }
    reference operator*() { return (*this)[0]; }
    const_reference operator*() const { return (*this)[0]; }
    bool operator== (const Column_iterator& i) const { return _same_col(i) && _row_offset==i._row_offset; }
    Column_iterator& operator++() { ++_row_offset; return *this; }
    Column_iterator operator++ (int) { Column_iterator tmp = *this; ++_row_offset; return tmp; }
    Column_iterator& operator--() { --_row_offset; return *this; }
    Column_iterator operator-- (int) { Column_iterator tmp = *this; --_row_offset; return tmp; }
    Column_iterator& operator+= (ptrdiff_t d) { _row_offset += d; return *this; }
    Column_iterator& operator-= (ptrdiff_t d) { _row_offset -= d; return *this; }
    friend Column_iterator operator+ (const Column_iterator& l, ptrdiff_t d) { Column_iterator res = l; res += d; return res; }
    friend Column_iterator operator- (const Column_iterator& l, ptrdiff_t d) { Column_iterator res = l; res -= d; return res; }
    ptrdiff_t operator- (const Column_iterator& i) const
    {
      if (!_same_col(i)) THROW Standard_exception ("Can't subtract array2d::Column_iterator's of different types");
      return _row_offset - i._row_offset;
    }
    Column_iterator begin() const { return *this; }
    Column_iterator end() const { return Column_iterator (_array, _col, _array->rows()); }
  };

  // row & column views
  class Row_view
  {
  private:
    Row_iterator _begin;
    Row_iterator _end;
  public:
    Row_view (array2d_type& a, int row)
      : _begin (a, row, 0), _end (a, row, a.columns())
    { }

    Row_view (array2d_type& a, int row, int begin_col, int end_col)
      : _begin (a, row, begin_col), _end (a, row, end_col)
    { }

    const Row_iterator& begin() const { return _begin; }
    const Row_iterator& end() const { return _end; }
    // override operator[] to do range checking for debug compilation
    reference operator[] (int col)
    {
#ifdef DART_DEBUG
      if (col < 0 || col >= _end - _begin) DART_DEBUG_ERROR ("Row_view overflow");
#endif /* DART_DEBUG */
      return _begin[col];
    }
    const_reference operator[] (int col) const
    {
#ifdef DART_DEBUG
      if (col < 0 || col >= _end - _begin) DART_DEBUG_ERROR ("Row_view overflow");
#endif /* DART_DEBUG */
      return _begin[col];
    }
  };

  class Column_view
  {
  private:
    Column_iterator _begin;
    Column_iterator _end;
  public:
    Column_view (array2d_type& a, int col)
      : _begin (a, col, 0), _end (a, col, a.columns())
    { }

    Column_view (array2d_type& a, int col, int begin_row, int end_row)
      : _begin (a, col, begin_row), _end (a, col, end_row)
    { }

    const Column_iterator& begin() const { return _begin; }
    const Column_iterator& end() const { return _end; }
    // override operator[] to do range checking for debug compilation
    reference operator[] (int col)
    {
#ifdef DART_DEBUG
      if (col < 0 || col >= _end - _begin) DART_DEBUG_ERROR ("Column_view overflow");
#endif /* DART_DEBUG */
      return _begin[col];
    }
    const_reference operator[] (int col) const
    {
#ifdef DART_DEBUG
      if (col < 0 || col >= _end - _begin) DART_DEBUG_ERROR ("Column_view overflow");
#endif /* DART_DEBUG */
      return _begin[col];
    }
  };

  // row & column accessors

  const Row_view row (int row) const { return Row_view (*this, row); }
  const Row_view row (int row, int begin_col, int end_col) const { return Row_view (*this, row, begin_col, end_col); }

  Row_view row (int row) { return Row_view (*this, row); }
  Row_view row (int row, int begin_col, int end_col) { return Row_view (*this, row, begin_col, end_col); }

  const Column_view column (int col) const { return Column_view (*this, col); }
  const Column_view column (int col, int begin_row, int end_row) const { return Column_view (*this, col, begin_row, end_row); }
  
  Column_view column (int col) { return Column_view (*this, col); }
  Column_view column (int col, int begin_row, int end_row) { return Column_view (*this, col, begin_row, end_row); }

  // misc get & set methods

  void set_row (int row, const vector<value_type>& v)
  {
    for (int x = 0; x < xsize(); x++) (*this)(x,row) = v[x];
  }

  void set_column (int col, const vector<value_type>& v)
  {
    for (int y = 0; y < ysize(); y++) (*this)(col,y) = v[y];
  }

  vector<vector<value_type> > all_columns() const
  {
    vector<vector<value_type> > v (ysize(), vector<value_type> (xsize(), (value_type) 0));
    for (int x = 0; x < xsize(); x++)
      for (int y = 0; y < ysize(); y++)
	v[x][y] = ((array2d_type&) *this)(x,y);
    return v;
  }

  array2d_type transpose() const
    {
      array2d_type tmp (ysize(), xsize());
      for (int x = 0; x < xsize(); x++)
	for (int y = 0; y < ysize(); y++)
	  tmp(y,x) = (*this)(x,y);
      return tmp;
    }

  // constructors
  array2d() : _xsize(0), _ysize(0), _data() {}

  array2d (const array2d_type& a) :
    _xsize(a._xsize), _ysize(a._ysize), _data(a._data) {}
  
  array2d (int xsize, int ysize, const value_type& t) :
    _xsize(xsize), _ysize(ysize), _data(xsize, ysize, t)
    {}
  
  array2d (int xsize, int ysize) :
    _xsize(xsize), _ysize(ysize), _data(xsize, ysize, value_type())
    {}

  // i/o methods
  void write_rowcol (ostream& out) const
  {
    for (int x = 0; x < xsize(); x++)
      for (int y = 0; y < ysize(); y++)
	out << ((array2d_type&) *this) (x, y) << (y < ysize()-1 ? " " : "\n");
  }
  void read_rowcol (istream& in)
  {
    for (int x = 0; x < xsize(); x++)
      for (int y = 0; y < ysize(); y++)
	in >> (*this) (x, y);
  }
};

template <class T, class xy_container_type>
ostream& operator<< (ostream& out, const array2d<T,xy_container_type>& a)
{
  for (int y = 0; y < a.ysize(); y++)
    for (int x = 0; x < a.xsize(); x++)
      out << ((array2d<T,xy_container_type>&) a) (x, y) << (x < a.xsize()-1 ? " " : "\n");
  return out;
}

template <class T, class xy_container_type>
istream& operator>> (istream& in, array2d<T,xy_container_type>& a)
{
  for (int y = 0; y < a.ysize(); y++)
    for (int x = 0; x < a.xsize(); x++)
      in >> a (x, y);
  return in;
}

#endif
