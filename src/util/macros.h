// miscellaneous useful macros

#ifndef MACROS_INCLUDED
#define MACROS_INCLUDED

#include <ctype.h>
#include <iostream>

// Import standard namespace.
using namespace std;

// DART version
#define DART_VERSION_STRING PACKAGE_VERSION   /* set by configure script */

// TRUE and FALSE
#define TRUE 1
#define FALSE 0

// conditional char and string macros for true and false
#define YES_OR_NO(X) ((X) ? "yes" : "no")
#define TRUE_OR_FALSE(X) ((X) ? "true" : "false")
#define Y_OR_N(X) ((X) ? 'y' : 'n')
#define T_OR_F(X) ((X) ? 't' : 'f')

// Container iterator macros. There are two kinds:
// (1) template iterator macros, using "typename";
// (2) nontemplate iterator macros, using typedefs without "typename".

// Prefixes used by the iterator macros
#define DART_ITER_TYPE _type_
#define DART_END_ITER  _end_
#define DART_TMP_CNTNR _cntnr_

// (1) template iterator macros, using "typename"
#define template_for_contents(ContainerType, Container, Iterator) for ( typename ContainerType ::iterator Iterator = ( Container ).begin(), DART_END_ITER ## Iterator = ( Container ).end(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define template_for_const_contents(ContainerType, Container, Iterator) for ( typename ContainerType ::const_iterator Iterator = ( Container ).begin(), DART_END_ITER ## Iterator = ( Container ).end(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define template_for_reverse_contents(ContainerType, Container, Iterator) for ( typename ContainerType ::reverse_iterator Iterator = ( Container ).rbegin(), DART_END_ITER ## Iterator = ( Container ).rend(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define template_for_const_reverse_contents(ContainerType, Container, Iterator) for ( typename ContainerType ::const_reverse_iterator Iterator = ( Container ).rbegin(), DART_END_ITER ## Iterator = ( Container ).rend(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )


// (2) nontemplate iterator macros, using typedefs without "typename"
#define for_contents(ContainerType, Container, Iterator) for ( ContainerType ::iterator Iterator = ( Container ).begin(), DART_END_ITER ## Iterator = ( Container ).end(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define for_const_contents(ContainerType, Container, Iterator) for ( ContainerType ::const_iterator Iterator = ( Container ).begin(), DART_END_ITER ## Iterator = ( Container ).end(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define for_tmp_contents(ContainerType, ContainerExpr, Iterator) ContainerType DART_TMP_CNTNR ## Iterator = ContainerExpr; for ( ContainerType ::iterator Iterator = DART_TMP_CNTNR ## Iterator .begin(), DART_END_ITER ## Iterator = DART_TMP_CNTNR ## Iterator .end(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define for_reverse_contents(ContainerType, Container, Iterator) for ( ContainerType ::reverse_iterator Iterator = ( Container ).rbegin(), DART_END_ITER ## Iterator = ( Container ).rend(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define for_const_reverse_contents(ContainerType, Container, Iterator) for ( ContainerType ::const_reverse_iterator Iterator = ( Container ).rbegin(), DART_END_ITER ## Iterator = ( Container ).rend(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

#define for_tmp_reverse_contents(ContainerType, ContainerExpr, Iterator) ContainerType DART_TMP_CNTNR ## Iterator = ContainerExpr; for ( ContainerType ::reverse_iterator Iterator = DART_TMP_CNTNR ## Iterator .rbegin(), DART_END_ITER ## Iterator = DART_TMP_CNTNR ## Iterator .rend(); !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

// for_iterator macro
#define for_iterator(IteratorType, Iterator, Begin, End) for ( IteratorType Iterator = Begin, DART_END_ITER ## Iterator = End; !(Iterator == DART_END_ITER ## Iterator); ++ Iterator )

// "begin, end" iterator macro
#define begin_end(Container) ( Container ).begin(), ( Container ).end()

// Directory separator
#define DIR_SEP_CHAR '/'


// Stream tricks -- rarely if ever used
struct Eat_white_type
{
  friend istream& operator>> (istream& is, Eat_white_type)
    {
      char c;
      while (is.get(c)) if (!isspace (c)) { is.putback(c); break; }
      return is;
    }
};

#define EAT_WHITE Eat_white_type()

struct Flush_type { friend ostream& operator<< (ostream& o, Flush_type) { return flush(o); } };

#define FLUSH Flush_type()

#ifdef DART_DEBUG
#define DART_DEBUG_MODE_STRING "debug"
#else
#define DART_DEBUG_MODE_STRING "release"
#endif

#endif
