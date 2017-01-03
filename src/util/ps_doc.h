#ifndef PS_DOC_INCLUDED
#define PS_DOC_INCLUDED

#include <stdio.h>
#include <iostream>
#include "util/sstring.h"

class PS_doc
{
 private:
  char*   _buf;
  sstring  _text;
  sstring  _header;
  sstring  _footer;

 public:
  PS_doc() :
    _buf (new char[1000]),
    _header ("%!PS-Adobe-3.0\n%%Pages: 1\n%%DocumentData: Clean7Bit\n%%PageOrder: Ascend\n%%Orientation: Portrait\n%%DocumentMedia: Default 612 792 0 () ()\n%%EndComments\n"),
    _footer ("%%Trailer\n%%bad%% end\n%%EOF\n")
    {}

  ~PS_doc() { delete[] _buf; }

  friend ostream& operator<<(ostream& o, const PS_doc& doc) { o << doc._header << doc._text << doc._footer; return o; }

  PS_doc& append (const char* txt) { _text.append(txt); return *this; }
  PS_doc& append () { return append (_buf); }

  sstring& header() { return _header; }
  sstring& footer() { return _footer; }

  PS_doc& scale (double xscale, double yscale)
    { sprintf (_buf, "%f %f scale\n", xscale, yscale); return append(); }

  PS_doc& translate (double dx, double dy)
    { sprintf (_buf, "%f %f translate\n", dx, dy); return append(); }

  PS_doc& rotate (double angle)
    { sprintf (_buf, "%f rotate\n", angle); return append(); }

  PS_doc& gsave () { return append ("gsave\n"); }
  PS_doc& grestore () { return append ("grestore\n"); }

  PS_doc& setgray (double g) { sprintf (_buf, "%f setgray\n", g); return append(_buf); }
  PS_doc& setlinewidth (double w) { sprintf (_buf, "%f setlinewidth\n", w); return append(_buf); }

  PS_doc& moveto (double x, double y)
    { sprintf (_buf, "%f %f moveto\n", x, y); return append(); }

  PS_doc& lineto (double x, double y)
    { sprintf (_buf, "%f %f lineto\n", x, y); return append(); }

  PS_doc& rlineto (double x, double y)
    { sprintf (_buf, "%f %f rlineto\n", x, y); return append(); }

  PS_doc& newpath () { return append ("newpath\n"); }
  PS_doc& closepath () { return append ("closepath\n"); }
  PS_doc& stroke () { return append ("stroke\n"); }
  PS_doc& fill () { return append ("fill\n"); }

  PS_doc& line (double x1, double y1, double x2, double y2) { moveto (x1, y1); lineto (x2, y2); return stroke(); }

  PS_doc& findfont (const char* font_name)
    { sprintf (_buf, "/%s findfont\n", font_name); return append(); }

  PS_doc& setfont () { return append("setfont\n"); }

  PS_doc& scalefont (double scale)
    { sprintf (_buf, "%f scalefont\n", scale); return append(); }

  PS_doc& show (const char* text)
    { sprintf (_buf, "(%s) show\n", text); return append(); }

  PS_doc& showpage () { return append ("showpage\n"); }

  PS_doc& cliprect (double x1, double y1, double x2, double y2)
    {
      sprintf (_buf, "%%clipRect\ninitclip\n%f %f moveto %f %f lineto %f %f lineto %f %f lineto %f %f lineto closepath eoclip newpath\n", x1, y1, x2, y1, x2, y2, x1, y2, x1, y1);
      return append (_buf);
    }

  PS_doc& setrgbcolor (double r, double g, double b)
    { sprintf (_buf, "%f %f %f setrgbcolor\n", r, g, b); return append (_buf); }

  PS_doc& rectangle (double x1, double y1, double x2, double y2)
    {
      return newpath().moveto(x1,y1).lineto(x2,y1).lineto(x2,y2).lineto(x1,y2).closepath();
    }
};

#endif
