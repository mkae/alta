;; Per-directory local variables for GNU Emacs 23 and later.

((python-mode
  ;; For Python, we want to display tabs as two columns, and we
  ;; basically want to use only tabs for indentation.
  . ((tab-width . 2)
     (indent-tabs-mode . t)
     (python-indent-offset . 2)))
 (c-mode
  . ((tab-width . 2)
     (indent-tabs-mode . t)
     (eval . (c-set-style "stroustrup"))))
 (c++-mode
  . ((tab-width . 2)
     (indent-tabs-mode . t)
     (eval . (c-set-style "stroustrup")))))
