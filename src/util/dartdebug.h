#ifndef DART_DEBUG_INCLUDED
#define DART_DEBUG_INCLUDED

#define DART_DEBUG
#define DART_DEBUG_ERROR(T) Dart_debug::runtime_error(T)

struct Dart_debug
{
  static void runtime_error (const char* text);
};

#endif
