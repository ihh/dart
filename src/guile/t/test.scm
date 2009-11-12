(define g (quote
	   ((grammar
	     (name JukesCantor)

	     ;; Transformation rules. These follow the pattern for a null model with rate matrix X.
	     (transform (from (S)) (to (X S*)))
	     (transform (from (S*)) (to (S)) (prob 1))
	     (transform (from (S*)) (to ()) (prob 1))
	     (update-rules 0)

	     ;; Rate matrix
	     (chain
	      ;; The state of this chain is a single nucleotide X.
	      (terminal (X))
	      ;; Treat initial probabilities and mutation rates as functions, not variables.
	      (update-policy parametric)

	      ;; initial probability distribution
	      (initial (state (a)) (prob .25))
	      (initial (state (c)) (prob .25))
	      (initial (state (g)) (prob .25))
	      (initial (state (t)) (prob .25))

	      ;; mutation rates
	      (mutate (from (a)) (to (c)) (rate .3333))
	      (mutate (from (a)) (to (g)) (rate .3333))
	      (mutate (from (a)) (to (t)) (rate .3333))
	      (mutate (from (c)) (to (a)) (rate .3333))
	      (mutate (from (c)) (to (g)) (rate .3333))
	      (mutate (from (c)) (to (t)) (rate .3333))
	      (mutate (from (g)) (to (a)) (rate .3333))
	      (mutate (from (g)) (to (c)) (rate .3333))
	      (mutate (from (g)) (to (t)) (rate .3333))
	      (mutate (from (t)) (to (a)) (rate .3333))
	      (mutate (from (t)) (to (c)) (rate .3333))
	      (mutate (from (t)) (to (g)) (rate .3333))))

	    ;; Standard DNA alphabet
	    (alphabet
	     (name DNA)
	     (token (a c g t))
	     (complement (t g c a))
	     (extend (to n) (from a) (from c) (from g) (from t))
	     (extend (to x) (from a) (from c) (from g) (from t))
	     (extend (to u) (from t))
	     (extend (to r) (from a) (from g))
	     (extend (to y) (from c) (from t))
	     (extend (to m) (from a) (from c))
	     (extend (to k) (from g) (from t))
	     (extend (to s) (from c) (from g))
	     (extend (to w) (from a) (from t))
	     (extend (to h) (from a) (from c) (from t))
	     (extend (to b) (from c) (from g) (from t))
	     (extend (to v) (from a) (from c) (from g))
	     (extend (to d) (from a) (from g) (from t))
	     (wildcard *)))))


(define s (stockholm-from-file "t/tiny.stk"))

(define w (xrate-train-grammar s g))  ;; should cause a warning
(define x (xrate-annotate-alignment s g))  ;; should cause a warning

(define t (xrate-estimate-tree s g))
(define u (xrate-train-grammar t g))
(define v (xrate-annotate-alignment t g))

(define ts (stockholm-tree s))  ;; should cause a warning
(define tt (stockholm-tree t))
(define ta (newick-ancestor-list tt))
(define tl (newick-leaf-list tt))
(define tb (newick-branch-list tt))

