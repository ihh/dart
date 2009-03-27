// Regexp.h - regular expression class based on Henry Spencer's regexp code.
// Adapted for Windows by Guy Gascoigne (see notes in Regexp.cc)
// Made unix-friendly by Ian Holmes ihh@lanl.gov April 30, 1999

#ifndef __REGEXP_H__
#define __REGEXP_H__

#include "util/sstring.h"

class regexp_internal;       //  internal interface to the regexp

class Regexp                 //  main regular expression class
{
public:
	enum { NSUBEXP = 10 };

	Regexp();
	Regexp( const char* exp, bool iCase = 0 );
	Regexp( const Regexp &r );
	~Regexp();
	const Regexp & operator=( const Regexp & r );

	bool Match( const char * s );
	int SubStrings() const;
	
	const sstring operator[]( unsigned int i ) const;
	int SubStart( unsigned int i ) const;
	int SubLength( unsigned int i ) const;

	sstring GetReplaceString( char* source ) const;

	sstring GetErrorString() const;
	bool CompiledOK() const;

#if defined( _RE_DEBUG )
	void Dump();
#endif
private:
	const char * str;	/* used to return substring offsets only */
	mutable sstring m_szError;
	regexp_internal * rc;

	void ClearErrorString() const;
	int safeIndex( unsigned int i ) const;

};

#endif

