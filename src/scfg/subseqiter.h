#ifndef SUBSEQITER_INCLUDED
#define SUBSEQITER_INCLUDED

#include "scfg/foldenv.h"

struct Subseq_iterator
{
  // data
  const Fold_envelope& xenv;
  const Fold_envelope& yenv;
  const Pair_envelope& pair_env;
  const int xsubseqs;
  const int ysubseqs;
  int xsubseq_idx;
  int ysubseq_idx;
  // constructors
  inline Subseq_iterator (const Fold_envelope& xenv, const Fold_envelope& yenv, const Pair_envelope& pair_env, int xsubseq_idx, int ysubseq_idx);
  inline Subseq_iterator (const Subseq_iterator& iter);
  // equality test operator
  inline int operator== (const Subseq_iterator& iter) const;
  inline int operator!= (const Subseq_iterator& iter) const;
  // accessors
  inline const Subseq& xsubseq() const;
  inline const Subseq& ysubseq() const;
protected:
  // protected inc/dec (called by ++ and -- operators in derived classes)
  inline void step_in();
  inline void step_out();
};

struct Inside_subseq_iterator : Subseq_iterator
{
  // constructors
  inline Inside_subseq_iterator (const Fold_envelope& xenv, const Fold_envelope& yenv, const Pair_envelope& pair_env, int xsubseq_idx, int ysubseq_idx);
  inline Inside_subseq_iterator (const Inside_subseq_iterator& iter);
  // inc/dec
  inline Inside_subseq_iterator& operator++();  // calls step_in()
  inline Inside_subseq_iterator operator++(int);
  inline Inside_subseq_iterator& operator--();  // calls step_out()
  inline Inside_subseq_iterator operator--(int);
};

struct Outside_subseq_iterator : Subseq_iterator
{
  // constructors
  inline Outside_subseq_iterator (const Fold_envelope& xenv, const Fold_envelope& yenv, const Pair_envelope& pair_env, int xsubseq_idx, int ysubseq_idx);
  inline Outside_subseq_iterator (const Outside_subseq_iterator& iter);
  // inc/dec (opposite direction to Inside_subseq_iterator)
  inline Outside_subseq_iterator& operator++();  // calls step_in()
  inline Outside_subseq_iterator operator++(int);
  inline Outside_subseq_iterator& operator--();  // calls step_out()
  inline Outside_subseq_iterator operator--(int);
};

// inline method defs

// Subseq_iterator
inline Subseq_iterator::Subseq_iterator (const Fold_envelope& xenv, const Fold_envelope& yenv, const Pair_envelope& pair_env, int xsubseq_idx, int ysubseq_idx)
  : xenv(xenv), yenv(yenv), pair_env(pair_env),
    xsubseqs(xenv.subseq.size()), ysubseqs(yenv.subseq.size()),
    xsubseq_idx(xsubseq_idx), ysubseq_idx(ysubseq_idx)
{ }

inline Subseq_iterator::Subseq_iterator (const Subseq_iterator& iter)
  : xenv(iter.xenv), yenv(iter.yenv), pair_env(iter.pair_env),
    xsubseqs(iter.xsubseqs), ysubseqs(iter.ysubseqs),
    xsubseq_idx(iter.xsubseq_idx), ysubseq_idx(iter.ysubseq_idx)
{ }

inline int Subseq_iterator::operator== (const Subseq_iterator& iter) const
{ return xsubseq_idx == iter.xsubseq_idx && ysubseq_idx == iter.ysubseq_idx; }

inline int Subseq_iterator::operator!= (const Subseq_iterator& iter) const
{ return xsubseq_idx != iter.xsubseq_idx || ysubseq_idx != iter.ysubseq_idx; }

inline void Subseq_iterator::step_in()
{
  do
    {
      if (++ysubseq_idx >= ysubseqs)
	{
	  ysubseq_idx = 0;
	  if (++xsubseq_idx >= xsubseqs) break;
	}
    }
  while (!Subseq_coords::allowed_by_pair_env (pair_env, xsubseq(), ysubseq()));
}

inline void Subseq_iterator::step_out()
{
  do
    {
      if (--ysubseq_idx < 0)
	{
	  ysubseq_idx = ysubseqs - 1;
	  if (--xsubseq_idx < 0) break;
	}
    }
  while (!Subseq_coords::allowed_by_pair_env (pair_env, xsubseq(), ysubseq()));
}

inline const Subseq& Subseq_iterator::xsubseq() const
{ return xenv.subseq[xsubseq_idx]; }

inline const Subseq& Subseq_iterator::ysubseq() const
{ return yenv.subseq[ysubseq_idx]; }

// Inside_subseq_iterator
inline Inside_subseq_iterator::Inside_subseq_iterator (const Fold_envelope& xenv, const Fold_envelope& yenv, const Pair_envelope& pair_env, int xsubseq_idx, int ysubseq_idx)
  : Subseq_iterator (xenv, yenv, pair_env, xsubseq_idx, ysubseq_idx)
{ }

inline Inside_subseq_iterator::Inside_subseq_iterator (const Inside_subseq_iterator& iter)
  : Subseq_iterator (iter)
{ }

inline Inside_subseq_iterator& Inside_subseq_iterator::operator++()
{ step_in(); return *this; }

inline Inside_subseq_iterator& Inside_subseq_iterator::operator--()
{ step_out(); return *this; }

inline Inside_subseq_iterator Inside_subseq_iterator::operator++(int)
{ Inside_subseq_iterator tmp (*this); step_in(); return tmp; }

inline Inside_subseq_iterator Inside_subseq_iterator::operator--(int)
{ Inside_subseq_iterator tmp (*this); step_out(); return tmp; }

// Outside_subseq_iterator
inline Outside_subseq_iterator::Outside_subseq_iterator (const Fold_envelope& xenv, const Fold_envelope& yenv, const Pair_envelope& pair_env, int xsubseq_idx, int ysubseq_idx)
  : Subseq_iterator (xenv, yenv, pair_env, xsubseq_idx, ysubseq_idx)
{ }

inline Outside_subseq_iterator::Outside_subseq_iterator (const Outside_subseq_iterator& iter)
  : Subseq_iterator (iter)
{ }

inline Outside_subseq_iterator& Outside_subseq_iterator::operator++()
{ step_out(); return *this; }

inline Outside_subseq_iterator& Outside_subseq_iterator::operator--()
{ step_in(); return *this; }

inline Outside_subseq_iterator Outside_subseq_iterator::operator++(int)
{ Outside_subseq_iterator tmp (*this); step_out(); return tmp; }

inline Outside_subseq_iterator Outside_subseq_iterator::operator--(int)
{ Outside_subseq_iterator tmp (*this); step_in(); return tmp; }

#endif /* SUBSEQITER_INCLUDED */
