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

typedef struct ss_elmt_str{
   char c; /*character either '<', '.', or '>'*/
   int index; /*element's index in stockholm ss_cons*/
   struct ss_elmt_str* under; /*pointer to last item put on stack*/
} ss_elmt;

typedef struct ss_stack_str{
   ss_elmt* top;
} ss_stack;

/*Functions to manipulate SS queues*/
ss_elmt* new_ss_elmt(char c, int i);
ss_stack* new_ss_stack();
ss_elmt* pop(ss_stack* s);
void push(ss_stack* s, ss_elmt* e);
int isempty(ss_stack* s);

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
