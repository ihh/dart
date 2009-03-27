#include <iostream>

using namespace std;

// had some very strange behaviour here --
// compiling via 'make', from emacs or command-line, gave an executable that does not catch the exception;
// however, typing the 'gcc' commands directly from the command-line, exactly as they are invoked by 'make',
// gave an executable that catches the exception just fine -- ihh, 26 Sept 2003

// It turned out that 'make' was using a different version of 'gcc'.
// 'gcc' was aliased to /usr/bin/gcc in ~/.cshrc, but make wasn't seeing this.
// ho hum.

// main
int main (int argc, char** argv)
{
  cerr << "(testing exception-catching)\n";
  try
    {
      throw 100;
    }
  catch (...)
    {
      cerr << "(caught exception: that's good)\n";
    }

  cout << "ok\n";
  return 0;
}
