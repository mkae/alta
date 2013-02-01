TEMPLATE = subdirs
SUBDIRS  = rational_function        \
           rational_data            \ 
			  rational_fitter_cgal     \
			  rational_fitter_quadprog \
			  rational_fitter_eigen

rational_fitter_cgal.depends     = rational_function rational_data
rational_fitter_quadprog.depends = rational_function rational_data
rational_fitter_eigen.depends    = rational_function rational_data
