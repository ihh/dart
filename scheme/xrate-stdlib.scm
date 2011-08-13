;; xrate standard library
;; Ian Holmes, August 13, 2011

;; Define fold, the most basic map-reduce functional construct that we need
;; (fold F (L1 L2 ... LN)) = (F L1 (F L2 ... (F LN '())))
;; This should be defined in SRFI-1, but that may not be around, so we use cond-expand
(cond-expand (srfi-1)
	     (else
	      (define fold
		(lambda (func list)
		  (if (equal? list '())
		      list
		      (func (car list) (fold func (cdr list))))))))

;; define ourselves a grep function
(define grep
  (lambda (predicate? list)
    (fold (lambda (item list) (if (predicate? item) (cons item list) list))
	  list)))

;; define a grep-not function
(define grep-not
  (lambda (predicate? list)
    (grep (lambda (item) (not (predicate? item))) list)))

;; define ourselves a Perl-style (spliced) map function
(define perl-map
  (lambda (func list)
    (fold (lambda (item fold-list) (append (func item) fold-list)) list)))

;; transition lookups
(define base-transition-table
  '((a g) (c t) (g a) (t c)))

;; genetic code lookups
(define codon-translation-table
  '(((a a a) K)
    ((a a c) N)
    ((a a g) K)
    ((a a t) N)
    ((a c a) T)
    ((a c c) T)
    ((a c g) T)
    ((a c t) T)
    ((a g a) R)
    ((a g c) S)
    ((a g g) R)
    ((a g t) S)
    ((a t a) I)
    ((a t c) I)
    ((a t g) M)
    ((a t t) I)
    ((c a a) Q)
    ((c a c) H)
    ((c a g) Q)
    ((c a t) H)
    ((c c a) P)
    ((c c c) P)
    ((c c g) P)
    ((c c t) P)
    ((c g a) R)
    ((c g c) R)
    ((c g g) R)
    ((c g t) R)
    ((c t a) L)
    ((c t c) L)
    ((c t g) L)
    ((c t t) L)
    ((g a a) E)
    ((g a c) D)
    ((g a g) E)
    ((g a t) D)
    ((g c a) A)
    ((g c c) A)
    ((g c g) A)
    ((g c t) A)
    ((g g a) G)
    ((g g c) G)
    ((g g g) G)
    ((g g t) G)
    ((g t a) V)
    ((g t c) V)
    ((g t g) V)
    ((g t t) V)
    ((t a a) !)
    ((t a c) Y)
    ((t a g) !)
    ((t a t) Y)
    ((t c a) S)
    ((t c c) S)
    ((t c g) S)
    ((t c t) S)
    ((t g a) !)
    ((t g c) C)
    ((t g g) W)
    ((t g t) C)
    ((t t a) M)
    ((t t c) M)
    ((t t g) M)
    ((t t t) M)))

;; (translation '(x y z)) returns the translation of codon xyz
(define translation
  (lambda (codon)
    (let ((codon-trans (assoc codon codon-translation-table)))
      (if (pair? codon-trans)
	  (cadr codon-trans)
	  #f))))

;; (transition x) returns the transition mutant of base x
(define transition
  (lambda (base)
    (let ((base-trans (assoc base base-transition-table)))
      (if (pair? base-trans)
	  (cadr base-trans)
	  #f))))

;; (one-mismatch? valid-mismatch? x y)
;; returns true if x and y differ at exactly one position, satisfying (valid-mismatch? x_value y_value)
(define one-mismatch?
  (lambda (valid-mismatch? x y)
    (if (list? x)
	(if (equal? x '())
	    #f
	    (or
	     (and (equal? (car x) (car y)) (one-mismatch? valid-mismatch? (cdr x) (cdr y)))
	     (and (one-mismatch? valid-mismatch? (car x) (car y)) (equal? (cdr x) (cdr y)))))
	(valid-mismatch? x y))))

(define any-one-mismatch?
  (lambda (x y)
    (one-mismatch? (lambda (a b) (not (equal? a b))) x y)))

;; (transition? x y) returns true if x & y (which may be lists or atoms) differ by a single-base transition
(define transition?
  (lambda (x y)
    (one-mismatch? (lambda (xnuc ynuc) (equal? (transition xnuc) ynuc)) x y)))

;; (transversion? x y) returns true if x & y (which may be lists or atoms) differ by a single-base transversion
(define transversion?
  (lambda (x y)
    (one-mismatch? (lambda (xnuc ynuc) (and (not (equal? xnuc ynuc)) (not (equal? (transition xnuc) ynuc)))) x y)))

;; (synonymous? (a b c) (d e f)) returns true if abc and def are synonymous codons
(define synonymous?
  (lambda (codon1 codon2)
    (equal? (translation codon1) (translation codon2))))


;; (stop-codon? (x y z)) returns true if xyz is a stop codon
(define stop-codon?
  (lambda (codon)
    (equal? (translation codon) '!)))


;; all-codons returns a list of all codons
(define all-codons
  (lambda ()
    (map (lambda (arg) (car arg))
	 codon-translation-table)))

;; all-stop-codons returns a list of all stop codons
(define all-stop-codons
  (lambda ()
    (grep stop-codon? (all-codons))))


;; all-sense-codons returns a list of all sense (non-stop) codons
(define all-sense-codons
  (lambda ()
    (grep-not stop-codon? (all-codons))))


;; helper function to test for zero/empty values
(define iszero? (lambda (x) (or (equal? x 0) (equal? x #f) (equal? x '()))))

;; function to build an xrate chain
;; (initfun x) should return initial probability of x
;; (mutfun x y) should return x->y mutation rate
(define xrate-chain
  (lambda (terminals states initfun mutfun)
    (quasiquote (chain
		 (update-policy parametric)
		 (terminal (unquote terminals))
		 (unquote-splicing (perl-map (lambda (state)
					       (let ((prob (initfun state)))
						 (if (iszero? prob)
						     '()
						     (quasiquote ((initial (state (unquote state)) (prob (unquote prob))))))))
					     states))
		 (unquote-splicing (perl-map (lambda (src)
					       (perl-map (lambda (dest)
							   (let ((rate (mutfun src dest)))
							     (if (iszero? rate)
								 '()
								 (quasiquote ((mutate (from (unquote src)) (to (unquote dest)) (rate (unquote rate)))))))) states)) states))))))

;; Codon chain
(define xrate-codon-chain
  (lambda (terminals initfun mutfun)
    (xrate-chain terminals (all-sense-codons) initfun mutfun)))

;; Simple parametric codon chain (single-substitution, transition/transversion, synonymous/nonsynonymous, c.f. Nielsen-Yang)
(define xrate-simple-codon-chain
  (lambda (terminals initfun transition-rate transversion-rate syn-rate nonsyn-rate)
    (xrate-codon-chain
     terminals
     initfun
     (lambda (src dest)
       (if (any-one-mismatch? src dest)
	   (let ((syn-factor (if (synonymous? src dest) syn-rate nonsyn-rate))
		 (trans-factor (if (transition? src dest) transition-rate transversion-rate))
		 (dest-factor (initfun dest)))
;; Uncomment to show multiplication signs explicitly
	     (quasiquote ((unquote trans-factor)
;;			  *
			  (unquote syn-factor)
;;			  *
			  (unquote dest-factor))))
	   #f)))))

;; codon-indexed parameters
(define xrate-codon-param
  (lambda (prefix codon)
    (string->symbol (string-concatenate (map symbol->string (cons prefix codon))))))

;; Nielsen-Yang model (slightly generalized)
(define xrate-NY-codon-chain
  (lambda (terminals eqm-prefix transition-rate transversion-rate syn-rate nonsyn-rate)
    (xrate-simple-codon-chain
     terminals
     (lambda (codon) (xrate-codon-param eqm-prefix codon))
     transition-rate transversion-rate syn-rate nonsyn-rate)))

;; probability parameter declarations for Nielsen-Yang model
(define xrate-NY-prob-params
  (lambda (eqm-prefix)
    (let ((num-codons (length (all-sense-codons))))
      (quasiquote (pgroup (unquote (map (lambda (codon)
					  (quasiquote ((unquote (xrate-codon-param eqm-prefix codon)) (unquote (/ 1 num-codons)))))
					(all-sense-codons))))))))

;; rate parameter declarations for Nielsen-Yang model
(define xrate-NY-rate-params
  (lambda rate-param-list
    (quasiquote (rate (unquote-splicing (map (lambda (param)
					       (quasiquote ((unquote param) 1)))
					     rate-param-list))))))

;; DNA alphabet
(define xrate-dna-alphabet '(alphabet (name DNA) (token (a c g t)) (complement (t g c a))))

;; Nielsen-Yang model as an xrate phylo-grammar
;; We leave this as a function so that the genetic code can be changed
(define xrate-NY-grammar
  (lambda ()
    (quasiquote
     (grammar
      (name Nielsen-Yang)
      (parametric)
      (transform (from (START)) (to ()) (prob 1))
      (transform (from (START)) (to (EMIT)) (prob 1))
      (transform (from (EMIT)) (to (NY1 NY2 NY3 EMIT*)) (prob 1))
      (transform (from (EMIT*)) (to (START)) (prob 1))
      (unquote (xrate-NY-prob-params 'pi_))
      (unquote (xrate-NY-rate-params 'kappa 'omega))
      (unquote (xrate-NY-codon-chain '(NY1 NY2 NY3) 'pi_ 'kappa 1 'omega 1))))))

;; Nielsen-Yang model as an xrate (alphabet, phylo-grammar) pair
;; As with xrate-NY-grammar, this is left as a function so that the genetic code can be changed
(define xrate-NY-alphgram
  (lambda ()
    (list xrate-dna-alphabet (xrate-NY-grammar))))
