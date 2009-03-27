#ifndef TELEGRAPH_ACCUMULATOR_INCLUDED
#define TELEGRAPH_ACCUMULATOR_INCLUDED

/* Vectors: extensible arrays with functions for easy resizing & reallocation.
   Entries -1 & -2 of the array are reserved for length & storage size.
*/

#define tg_vector(type)  type*
typedef tg_vector(void*) tg_vec;

/* tg_alloc: lvalue for allocated storage of a tg_vector (entry -2)
   tg_size: lvalue for current size of a tg_vector (entry -1)
 */
#define tg_alloc(vec) (((unsigned int*) vec)[-2])
#define tg_size(vec)  (((unsigned int*) vec)[-1])

/* Functions for vectors.
   Many of these return the new (possibly changed) tg_vec.
 */
tg_vec tg_vec_new (int initAlloc);            /* constructor */
tg_vec tg_vec_init (int initSize, void* val); /* constructor/initialiser */
void tg_vec_delete (tg_vec* vec);             /* destructor, sets vec = null */
int tg_vec_append (tg_vec* vec, void* value); /* returns new element index */
void tg_vec_reserve (tg_vec* vec, int alloc); /* ensures storage >= alloc */
void tg_vec_clear (tg_vec vec);               /* sets size to 0 */
int tg_vecs_equal (tg_vec vec1, tg_vec vec2); /* compares sizes & elements */

/* Basic vector types */
typedef tg_vector(unsigned int) tg_unsigned_vec;
typedef tg_vector(const char*)  tg_string_vec;
typedef tg_vector(unsigned int) tg_char_vec; 

#endif /* TELEGRAPH_ACCUMULATOR_INCLUDED */

