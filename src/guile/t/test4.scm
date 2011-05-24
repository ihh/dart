;; Running 'darts -s <this file>' and 'terminatrix -l -r ../../ontology/t/test6.tsm' should generate identical output

;; uncomment these (dart-log...) lines for more debugging output...
;; (dart-log "TERMINATRIX")
;; (dart-log "GUILE")
;; (dart-log "5")

(define
  my_tree_db
  (list (list "my_tree_name" (newick-from-string "(A:1,B:1);"))))

(define
  (my_knowledge family gene state-tuple)
  (if
   (equal? gene "A")
   (equal? (car state-tuple) "term1")
   (equal? (car state-tuple) "term2")))

(define
  my_model
  (quote
   (model
    (chain
     (update-policy parametric)
     (terminal (TERM))

     ;; initial probability distribution
     (initial (state (term1)) (prob p1))
     (initial (state (term2)) (prob p2))
     (initial (state (term3)) (prob p3))
     (initial (state (term4)) (prob p4))

     ;; mutation rates
     (mutate (from (term1)) (to (term2)) (rate (r p2)))
     (mutate (from (term1)) (to (term3)) (rate (r p3)))
     (mutate (from (term1)) (to (term4)) (rate (r p4)))
     (mutate (from (term2)) (to (term1)) (rate (r p1)))
     (mutate (from (term2)) (to (term3)) (rate (r p3)))
     (mutate (from (term2)) (to (term4)) (rate (r p4)))
     (mutate (from (term3)) (to (term1)) (rate (r p1)))
     (mutate (from (term3)) (to (term2)) (rate (r p2)))
     (mutate (from (term3)) (to (term4)) (rate (r p4)))
     (mutate (from (term4)) (to (term1)) (rate (r p1)))
     (mutate (from (term4)) (to (term2)) (rate (r p2)))
     (mutate (from (term4)) (to (term3)) (rate (r p3))))

    (alphabet
     (name GO)
     (token (term1 term2 term3 term4))))))

(define
  my_pgroups
  (quote
   (pgroup
    ((p1 .25) (p2 .25) (p3 .25) (p4 .25)))))

(define
  my_rates
  (quote
   (rate
    (r 1.33))))

(define
  results
  (terminatrix-learn
   100
   (quasiquote
    (terminatrix
     (knowledge-scheme my_knowledge)
     (tree-db-scheme my_tree_db)
     (model-scheme my_model)
     (params (unquote-splicing (list my_pgroups my_rates)))))))

;; Mimic terminatrix output
(display ";; EM results\n")
(display results)
(display "\n")
