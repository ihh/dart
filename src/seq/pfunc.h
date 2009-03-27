#ifndef PFUNC_INCLUDED
#define PFUNC_INCLUDED

#include "seq/pvar.h"

// A PFunc is a general class for representing simple probability expressions, eg "(A+(B*C))/(D+E...)"

// PFuncs are stored as a sequence of operators in reverse Polish notation,
//  corresponding to simple algebraic expressions involving stack-pushes of constants and PVar's
//  and the "probabilistic" stack-reduction operations: '+' '*' and '/'.

// (Don't ask me about Borel spaces. I don't wanna know.)

class PFunc
{
public:
  // operator types
  // NB INERT_PVAR is the same as PVAR, except its derivative w.r.t. the corresponding PVar is zero
  enum Opcode { CONST, PVAR, INERT_PVAR, ADD, MUL, DIV, POW };

private:
  // reverse Polish expression
  vector<Opcode> opcode;    // the opcodes
  vector<Score>  const_sc;  // the scores (for CONST nodes)
  vector<PVar>   pvar;      // the PVar's (for PVAR and INERT_PVAR nodes)

  // Structure of internal pointers and node scores for the reverse Polish
  //  stack.  Used by eval_sc & eval_dlogf_dlogp
  struct Stack_index
  {
    const PFunc& f;              // the PFunc
    vector<int>    l_operand_idx;  // for ADD/MUL/DIV/POW nodes, the index of the left operand node
    vector<int>    r_operand_idx;  // for ADD/MUL/DIV nodes, the index of the right operand node
    vector<int>    cp_idx;         // for CONST/POW/PVAR/INERT_PVAR nodes, the index into the Score or PVar vectors, respectively
    vector<Score>  node_sc;        // evaluated Score's at each node

    Stack_index (const PFunc& f);
    void eval_node_sc (const PScores& var_sc);
  };
  friend struct Stack_index;

public:
  // main interface

  // method to test if this func is null
  bool is_null() const;

  // method to test if this func is obviously zero (NB not the same as null!)
  bool is_zero() const;

  // method to test if this func is obviously one
  bool is_one() const;

  // method to return all dependent variables
  set<PVar> vars() const;

  // method to return all dependent groups
  set<PGroup> groups() const;

  // evaluate Prob2Score(F)  where F is this function
  Score eval_sc (const PScores& var_sc) const;

  // evaluate d(log(F)) / d(log(pvar))
  Prob  eval_dlogf_dlogp (const PScores& var_sc,
			  const PVar& pvar) const;

  // incremement PVar counts
  void  inc_var_counts (PCounts& var_counts,
			const PScores& var_sc,
			Prob func_count,
			Prob func_wait = 0.) const;
  
  // builder methods
  // constructors
  PFunc();   // returns a null function
  PFunc (Prob constant);
  PFunc (PVar pvar);
  PFunc (Opcode op, PFunc fl, PFunc fr);
  PFunc (PFunc f, Prob exponent);

  // method to transform all PVAR opcodes into INERT_PVAR
  void make_inert();

  // mutating operators
  PFunc& operator= (Prob r)   { return *this = PFunc(r); }
  PFunc& operator= (PVar r)   { return *this = PFunc(r); }

  PFunc& operator+= (Prob r)  { return *this = PFunc (ADD, *this, PFunc(r)); }
  PFunc& operator+= (PVar r)  { return *this = PFunc (ADD, *this, PFunc(r)); }
  PFunc& operator+= (PFunc r) { return *this = PFunc (ADD, *this, PFunc(r)); }

  PFunc& operator*= (Prob r)  { return *this = PFunc (MUL, *this, PFunc(r)); }
  PFunc& operator*= (PVar r)  { return *this = PFunc (MUL, *this, PFunc(r)); }
  PFunc& operator*= (PFunc r) { return *this = PFunc (MUL, *this, PFunc(r)); }

  PFunc& operator/= (Prob r)  { return *this = PFunc (DIV, *this, PFunc(r)); }
  PFunc& operator/= (PVar r)  { return *this = PFunc (DIV, *this, PFunc(r)); }
  PFunc& operator/= (PFunc r) { return *this = PFunc (DIV, *this, PFunc(r)); }
  
  // non-mutating operators
  friend PFunc operator+ (Prob l,  PVar r)  { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (Prob l,  PFunc r) { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (PVar l,  Prob r)  { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (PVar l,  PVar r)  { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (PVar l,  PFunc r) { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (PFunc l, Prob r)  { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (PFunc l, PVar r)  { return PFunc (ADD, PFunc(l), PFunc(r)); }
  friend PFunc operator+ (PFunc l, PFunc r) { return PFunc (ADD, PFunc(l), PFunc(r)); }

  friend PFunc operator* (Prob l,  PVar r)  { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (Prob l,  PFunc r) { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (PVar l,  Prob r)  { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (PVar l,  PVar r)  { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (PVar l,  PFunc r) { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (PFunc l, Prob r)  { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (PFunc l, PVar r)  { return PFunc (MUL, PFunc(l), PFunc(r)); }
  friend PFunc operator* (PFunc l, PFunc r) { return PFunc (MUL, PFunc(l), PFunc(r)); }

  friend PFunc operator/ (Prob l,  PVar r)  { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (Prob l,  PFunc r) { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (PVar l,  Prob r)  { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (PVar l,  PVar r)  { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (PVar l,  PFunc r) { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (PFunc l, Prob r)  { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (PFunc l, PVar r)  { return PFunc (DIV, PFunc(l), PFunc(r)); }
  friend PFunc operator/ (PFunc l, PFunc r) { return PFunc (DIV, PFunc(l), PFunc(r)); }

  // display methods
  void show (ostream& o, const vector<vector<sstring> >* group_suffix = 0, bool no_group_prefix = false) const;
};

PFunc operator+ (Prob l,  PVar r);
PFunc operator+ (Prob l,  PFunc r);
PFunc operator+ (PVar l,  Prob r);
PFunc operator+ (PVar l,  PVar r);
PFunc operator+ (PVar l,  PFunc r);
PFunc operator+ (PFunc l, Prob r);
PFunc operator+ (PFunc l, PVar r);
PFunc operator+ (PFunc l, PFunc r);

PFunc operator* (Prob l,  PVar r);
PFunc operator* (Prob l,  PFunc r);
PFunc operator* (PVar l,  Prob r);
PFunc operator* (PVar l,  PVar r);
PFunc operator* (PVar l,  PFunc r);
PFunc operator* (PFunc l, Prob r);
PFunc operator* (PFunc l, PVar r);
PFunc operator* (PFunc l, PFunc r);

PFunc operator/ (Prob l,  PVar r);
PFunc operator/ (Prob l,  PFunc r);
PFunc operator/ (PVar l,  Prob r);
PFunc operator/ (PVar l,  PVar r);
PFunc operator/ (PVar l,  PFunc r);
PFunc operator/ (PFunc l, Prob r);
PFunc operator/ (PFunc l, PVar r);
PFunc operator/ (PFunc l, PFunc r);

// Kronecker delta class for PFuncs

class Kronecker_func : public vector<PFunc>
{
public:
  Kronecker_func (int size, int val) : vector<PFunc> (size, PFunc(0.0)) { (*this)[val] = PFunc(1.0); }
};


#endif
