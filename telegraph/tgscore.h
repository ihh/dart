#ifndef TELEGRAPH_SCORE_INCLUDED
#define TELEGRAPH_SCORE_INCLUDED

/* 80-column formatting
.........1.........2.........3.........4.........5.........6.........7.........
1234567890123456789012345678901234567890123456789012345678901234567890123456789
*/

  /* Scores */
  typedef int tg_score;

#define tg_score2bit_ratio 1000.
#define tg_infinite_score  987654321
#define tg_undefined_score 998877666
#define tg_log2_infinity   (((double) tg_infinite_score) / tg_score2bit_ratio)
#define tg_score2prob(S) (pow(2,((double)(S))/tg_score2bit_ratio))
#define tg_prob2score(P) ((tg_score)(tg_score2bit_ratio*log((P))/log(2)))
#define tg_score2bits(S)  (((double)(S))/tg_score2bit_ratio)

typedef tg_vec tg_score_array;  /* elements have type (tg_score) */

#endif /* TELEGRAPH_SCORE_INCLUDED */
