#include "handel/recorder.h"

// suffix for auto-saved dotfiles
#define DOTFILE_SUFFIX ".dot"

// suffix for auto-saved S-expressions
#define SEXPR_SUFFIX ".sxpr"

ofstream* Numbered_file_opener::open_file()
{
  sstring filename;
  filename << directory << DIR_SEP_CHAR << prefix << ++n_files << suffix;
  CTAG(4,RECORDER) << "Saving file to '" << filename << "'\n";

  ofstream* out = new ofstream (filename.c_str());
  if (!*out)
    THROWEXPR ("Couldn't open file '" << filename << "'");

  return out;
}

Transducer_dotfile_recorder::Transducer_dotfile_recorder()
  : Numbered_file_opener()
{
  suffix = DOTFILE_SUFFIX;
}

void Transducer_dotfile_recorder::save (const Transducer<Score>& trans, const char* comment)
{
  if (!active())
    return;

  ofstream* out = open_file();

  if (comment)
    *out << "// " << comment << "\n";
  trans.print_dotfile (*out);

  out->close();
  delete out;
}

Composition_recorder::Composition_recorder()
  : Numbered_file_opener()
{
  suffix = SEXPR_SUFFIX;
}

void Composition_recorder::save (Handel_movement& move)
{
  if (active())
    {
      ofstream* out = open_file();
      move.dump_composition (*out);
      out->close();
      delete out;
    }
  else
    {
      sstring dummy;
      move.dump_composition (dummy);
    }
}
