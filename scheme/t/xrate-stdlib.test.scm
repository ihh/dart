(primitive-load-path "xrate-stdlib.scm")

;; Tests
(display "translation (1 2 3) = ")
(display (translation '(1 2 3)))
(newline)

(display "translation (t a t) = ")
(display (translation (quote (t a t))))
(newline)

(display "translation (t a g) = ")
(display (translation (quote (t a g))))
(newline)


(display "stop-codon? (a g t) = ")
(display (stop-codon? (quote (a g t))))
(newline)

(display "stop-codon? (t g a) = ")
(display (stop-codon? (quote (t g a))))
(newline)


(display "transition t = ")
(display (transition (quote t)))
(newline)


(display "transition? a c = ")
(display (transition? (quote a) (quote c)))
(newline)

(display "transition? a g = ")
(display (transition? (quote a) (quote g)))
(newline)


(display "transition? (a a a) (a a a) = ")
(display (transition? (quote (a a a)) (quote (a a a))))
(newline)

(display "transition? (a a a) (a c a) = ")
(display (transition? (quote (a a a)) (quote (a c a))))
(newline)

(display "transition? (a a a) (a a g) = ")
(display (transition? (quote (a a a)) (quote (a a g))))
(newline)

(display "transition? (a a a) (a c g) = ")
(display (transition? (quote (a a a)) (quote (a c g))))
(newline)


(display "transversion? a c = ")
(display (transversion? (quote a) (quote c)))
(newline)

(display "transversion? a g = ")
(display (transversion? (quote a) (quote g)))
(newline)


(display "transversion? (a a a) (a a a) = ")
(display (transversion? (quote (a a a)) (quote (a a a))))
(newline)

(display "transversion? (a a a) (a c a) = ")
(display (transversion? (quote (a a a)) (quote (a c a))))
(newline)

(display "transversion? (a a a) (a a g) = ")
(display (transversion? (quote (a a a)) (quote (a a g))))
(newline)

(display "transversion? (a a a) (a c g) = ")
(display (transversion? (quote (a a a)) (quote (a c g))))
(newline)


(display "synonymous? (a a a) (a a c) = ")
(display (synonymous? (quote (a a a)) (quote (a a c))))
(newline)


(display "synonymous? (a a a) (a a g) = ")
(display (synonymous? (quote (a a a)) (quote (a a g))))
(newline)


(display "all-codons = ")
(display (all-codons))
(newline)

(display "all-stop-codons = ")
(display (all-stop-codons))
(newline)

(display "all-sense-codons = ")
(display (all-sense-codons))
(newline)

(display "xrate-NY-alphgram = ")
(display (xrate-NY-alphgram))
(newline)

(display "Validating grammar ... ")
(xrate-validate-grammar (xrate-NY-alphgram))
(newline)

(display 'done)
(newline)
