#ifndef DART_UNIX_ENVIRONMENT_INCLUDED
#define DART_UNIX_ENVIRONMENT_INCLUDED

#include <stdlib.h>

#define DARTDIR_ENV_VAR "DARTDIR"

struct Dart_Unix
{
  static char* get_DARTDIR() { return getenv (DARTDIR_ENV_VAR); }
};

#endif /* DART_UNIX_ENVIRONMENT_INCLUDED */
