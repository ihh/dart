// my exception class

#ifndef DART_EXCEPTION_INCLUDED
#define DART_EXCEPTION_INCLUDED

#include "util/sstring.h"

struct Dart_exception
{
  sstring stack_trace, _what;

  Dart_exception();
  virtual ~Dart_exception();

  virtual const char* details() const;
  const char* what() const;
};

class Standard_exception : public Dart_exception
{
 protected:
  sstring msg;

 public:
  Standard_exception (const char* m);
  Standard_exception (const sstring& m);
  const char* details() const;
};

class Format_exception : public Dart_exception
{
 protected:
  sstring info;
  void setup_info (istream& in);

 public:
  Format_exception (istream& in);
  Format_exception (istream& in, const char* prefix);
  Format_exception (istream& in, const sstring& prefix);
  const char* details() const;
};

class String_exception : public Dart_exception
{
 protected:
  sstring msg;

 public:
  String_exception (const char* m, const char* arg) : Dart_exception(), msg(m) { msg.append(arg).append("\n"); }
  String_exception (const char* m, const sstring& arg) : Dart_exception(), msg(m) { msg.append(arg).append("\n"); }
  const char* details() const;
};

#endif
