#ifndef NEWTON_RAPHSON_INCLUDED
#define NEWTON_RAPHSON_INCLUDED

#include <vector>
#include <stack>
#include <set>
#include <map>
#include "util/sstring.h"

// Grad_function is a virtual class for evaluating a function $f({\bf x})$ and its first derivatives
// w.r.t. a parameter vector $\bf x$.
struct Grad_function
{
  virtual int dim() = 0;  // size of parameter vector
  virtual void eval (const vector<double>& x, double& f, vector<double>& grad_f) = 0;  // evaluate function & derivatives
  virtual ~Grad_function() { }
};

// Newton_Raphson is a virtual class for doing Newton-Raphson on a Grad_function
// with a vector of constraints, $\bf a$, such that $x_i > a_i \forall i$.
struct Newton_Raphson
{
  typedef map<int,double> Constraints;
  static vector<double> iterate (Grad_function& f, const vector<double>& xseed, const Constraints& xmin, double tol, int max_iter);
};

// A Diff is a differentiable function
//  that uses a reverse Polish stack
class Diff
{
public:
  enum Opcode { CONST, PARAM, ADD, MUL, LOG, EXP, POW,  // these are the actual operators
		NOP, NOPC, NOPP };  // these are no-op placeholders (the latter two consume a const/param, respectively)

  const vector<sstring>* param_name;
  sstring make_param_name (int i) const;

protected:
  // reverse Polish expression
  vector<Opcode> opcode;    // the opcodes
  vector<double> constant;  // the constants (for CONST and POW nodes)
  vector<int>    param;     // the parameter indices (for PARAM nodes)

  // Structure of internal stack pointers
  struct Stack_index
  {
    // referenced Diff
    Diff& D;
    
    // parser variables
    stack<int> idx_stack;
    int const_i;
    int param_i;

    // lookup tables
    vector<int> l_op_idx;   // for ADD/MUL/LOG/EXP/POW nodes, the index of the left operand node
    vector<int> r_op_idx;   // for ADD/MUL nodes, the index of the right operand node (which happens to always be the previous node)
    vector<int> start_idx;  // index of the first opcode in this subexpression
    vector<int> const_idx;  // index of the next constant[]
    vector<int> param_idx;  // index of the next param[]

    // constructor
    Stack_index (Diff& D, bool _update = 1);

    // test methods
    void assert_single_stack() const;

    // update methods
    void forget (int expr_idx = 0);  // forgets everything from expr_idx onwards; invalidates idx_stack
    void update (int new_size);
    void update();

    // edit methods
    void insert_node (int expr_idx, Diff& dest, int opcode_pos, int const_pos, int param_pos) const;
    void append_node (int expr_idx, Diff& dest) const;
    void delete_node (int expr_idx);  // invalidates idx_stack
    void swap_last_two_nodes();
  };
  friend struct Stack_index;

  void append (const Diff& diff);
  bool is_null() const;

public:
  // constructor
  Diff() : param_name(0) { }

  // dimension
  int dim() const;

  // parameters
  set<int> params() const;

  // eval methods
  double eval (const vector<double>& x) const;
  Diff derivative (int deriv_idx) const;

  // builder methods
  void clear();
  
  Diff& operator*= (const Diff& diff);
  Diff& operator+= (const Diff& diff);
  Diff& operator-= (const Diff& diff);
  Diff& operator/= (const Diff& diff);

  friend Diff operator* (const Diff& a, const Diff& b);
  friend Diff operator+ (const Diff& a, const Diff& b);
  friend Diff operator- (const Diff& a, const Diff& b);
  friend Diff operator/ (const Diff& a, const Diff& b);

  Diff log() const;
  Diff exp() const;
  Diff pow (double p) const;

  static Diff sum (vector<Diff>& series);

  // simplify method: applies some very crude optimisations to the tree
  void simplify();

  // display methods
  sstring to_string() const;
  sstring params_to_string (const vector<double>& x) const;
};

// trivial subclasses of Diff
class Diff_const : public Diff { public: Diff_const (double c); };
class Diff_param : public Diff { public: Diff_param (int i); };
class Diff_log   : public Diff { public: Diff_log (const Diff& diff); };
class Diff_exp   : public Diff { public: Diff_exp (const Diff& diff); };
class Diff_power : public Diff { public: Diff_power (const Diff& diff, double p); };

// should probably have a Grad_function wrapper for Diff here.

// A Laplacian is a Grad_function of the form (grad D)^2, where D is a Diff
struct Laplacian : Grad_function
{
  const Diff& D;  // F = (grad D)^2
  vector<Diff> dD_dx;  // dD_dx[i] = D.derivative(i)
  typedef map<int,Diff> Diff_map;
  vector<Diff_map> d2D_dxdy;  // d2D_dxdy[i][j] = D.derivative(i).derivative(j), j <= i
  Laplacian (const Diff& D);
  int dim();
  void eval (const vector<double>& x, double& f, vector<double>& grad_f);
  void eval_cmp (const vector<double>& x, double& f, vector<double>& grad_f, double& cmp);  // sets cmp = D
};

#endif /* NEWTON_RAPHSON_INCLUDED */
