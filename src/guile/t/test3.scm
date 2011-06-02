;; test for ontology term evolution
(dart-log "TERMX")
(dart-log "GUILE")
(dart-log "5")

(define
  my_tree_db
  (list (list "my_tree_name" (newick-from-string "(A:1,B:1);"))))

(define
  (my_knowledge family gene state-tuple)
  (if
   (and
    (equal? gene "A")
    (not (equal? (car state-tuple) "term1")))
   #f
   #t))

(define
  my_model
  (quote
   (model
    (chain
     (update-policy parametric)
     (terminal (TERM))

     ;; initial probability distribution
     (initial (state (term1)) (prob .25))
     (initial (state (term2)) (prob .25))
     (initial (state (term3)) (prob .25))
     (initial (state (term4)) (prob .25))

     ;; mutation rates
     (mutate (from (term1)) (to (term2)) (rate .33))
     (mutate (from (term1)) (to (term3)) (rate .33))
     (mutate (from (term1)) (to (term4)) (rate .33))
     (mutate (from (term2)) (to (term1)) (rate .33))
     (mutate (from (term2)) (to (term3)) (rate .33))
     (mutate (from (term2)) (to (term4)) (rate .33))
     (mutate (from (term3)) (to (term1)) (rate .33))
     (mutate (from (term3)) (to (term2)) (rate .33))
     (mutate (from (term3)) (to (term4)) (rate .33))
     (mutate (from (term4)) (to (term1)) (rate .33))
     (mutate (from (term4)) (to (term2)) (rate .33))
     (mutate (from (term4)) (to (term3)) (rate .33)))

    (alphabet
     (name GO)
     (token (term1 term2 term3 term4))))))

(write
 (termx-evidence
  (quote
   (termx
    (knowledge-scheme my_knowledge)
    (tree-db-scheme my_tree_db)
    (model-scheme my_model)))))

(write
 (termx-predict
  (quote
   (termx
    (knowledge-scheme my_knowledge)
    (tree-db-scheme my_tree_db)
    (model-scheme my_model)))))
