#include "telegraph/tgvar.h"
#include "seq/biosequence.h"
#include "util/Regexp.h"

int main (int argc, char** argv)
{
  try {
    // first test the troublesome Regexp
    /* Commented out because this definitely causes a crash with the Henry Spencer Regexp library
       Regexp re1 = "^[ \t\n]*([ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_][ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789]*)[ \t\n]*=[ \t\n]*(.+)$";
       Regexp re2 = "^[ \t\n]*([ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_][ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789]*)";
       Regexp re2_5 = "([ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_][ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789]*)";
       Regexp re2_6 = "[ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_][ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789]*";
       Regexp re3 = "^[ \t\n]*([ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_])";
       Regexp re4 = "^[ \t\n]*";
       const char* cstr = "NullEmitA = 2**-2.672";

       if (re4.Match (cstr)) cerr << "(matched Regexp re4)\n"; else { cout << "not ok\n"; exit(1); }
       if (re3.Match (cstr)) cerr << "(matched Regexp re3)\n"; else { cout << "not ok\n"; exit(1); }
       if (re2_6.Match (cstr)) cerr << "(matched Regexp re2_6)\n"; else { cout << "not ok\n"; exit(1); }
       if (re2_5.Match (cstr)) cerr << "(matched Regexp re2_5)\n"; else { cout << "not ok\n"; exit(1); }
       if (re2.Match (cstr)) cerr << "(matched Regexp re2)\n"; else { cout << "not ok\n"; exit(1); }
       if (re1.Match (cstr)) cerr << "(matched Regexp re1)\n"; else { cout << "not ok\n"; exit(1); }
    */

    // now test the Telegraph_PScores_adaptor
    const Alphabet& alphabet = RNA_alphabet;
    PScores pscore;
    const Alphabet_group null_emit = pscore.new_alphabet_group (alphabet, "NullEmit");
    Telegraph_PScores_adaptor tgio (pscore);
    tgio.read ("TGVAR-TEST");
  }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << "ok\n";
}
