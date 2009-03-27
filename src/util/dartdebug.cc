#include <iostream>
#include "util/dartdebug.h"
#include "util/logfile.h"

void Dart_debug::runtime_error (const char* text)
{
  THROW Standard_exception (text);
}
