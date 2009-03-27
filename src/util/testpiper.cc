#include "util/piper.h"
#include "util/sstring.h"
#include "util/logfile.h"

using namespace std;

struct Piper_test : Piper
{
  Piper_test (int max_fork) : Piper (max_fork) { }
  void test()
  {
    open_pipes();
    for (int proc = 0; proc < max_fork; ++proc)
      if (fork_child(proc))
	{
	  // child
	  const int f = write_fd (proc);
	  cerr << "[emitting HelloWorld#" << proc << " on fd#" << f << "]\n";
	  ostream_cfile pipe (f);
	  pipe.out << "Hello\n";
	  pipe.out << "world\n";
	  pipe.out << proc << "\n";
	  pipe.out.flush();
	  pipe.close();
	  exit(0);
	}
      else
	sleep (1);  // parent; pause here
    for (int proc = 0; proc < max_fork; ++proc)
      {
	const int f = read_fd (proc);
	cerr << "[waiting for HelloWorld#" << proc << " on fd#" << f << "]\n";
	istream_cfile pipe (f);
	sstring s;
	s.getline (pipe.in);
	if (s != sstring ("Hello")) THROWEXPR ("Expected 'Hello'");
	s.getline (pipe.in);
	if (s != sstring ("World")) THROWEXPR ("Expected 'World'");
	int p;
	pipe.in >> p;
	if (p != proc) THROWEXPR ("Expected " << p);
      }
  }
};

int main (int argc, char** argv)
{
  Piper_test piper_test (2);
  try
    {
      piper_test.test();
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }
  cout << "ok\n";
  return 0;
}
