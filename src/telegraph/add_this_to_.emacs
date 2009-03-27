;; flex/bison
;;
(load "make-regexp.elc")
(load "flex-mode.el")
(load "bison-mode.el")
 
(setq auto-mode-alist (cons '("\\.l\\'" . flex-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.y\\'" . bison-mode) auto-mode-alist))                                                                                                            
