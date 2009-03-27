#ifndef PAM_INCLUDED
#define PAM_INCLUDED

#include "tree/diagonalised_matrix_factory.h"

struct PAM_factory : Diagonalised_matrix_factory
{
  PAM_factory();
  const Alphabet& alphabet();
};

#endif
