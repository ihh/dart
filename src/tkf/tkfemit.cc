#define __PUT_STATIC_DATA_MEMBERS_HERE

#include "tkf/tkfsequence.h"
#include "tkf/tkfopts.h"
#include "util/rnd.h"
#include "util/logfile.h"

int main(int argc, char* argv[])
{
  TKF_opts opts (argc, argv, 0);
  opts.short_description = "generate sequences related by a tree using the TKF model";
  opts.syntax = "[options] <treefile>";

  opts.newline();
  opts.print ("Sequence generation options\n");
  opts.print ("---------------------------\n");

  bool use_expected_length;
  opts.add ("expectedlen", use_expected_length = 0, "fix root sequence length at expected value");

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 1) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // display initial logging messages

      Score_fns::describe_scoring_scheme (CLOG(8));

      // do the work
      
      TKF_params params = opts.params();
      TKF_emit_align tkf (params);

      ifstream tree_file (opts.args[0].c_str());
      if (!tree_file) THROW String_exception ("Tree file not found: ", opts.args[0].c_str());

      CLOG(7) << "Reading tree\n";
      tkf.read_PHYLIP (tree_file);
      
      tkf.emit (use_expected_length);
      tkf.align.write_MUL (cout, params.submat_factory.alphabet());
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
