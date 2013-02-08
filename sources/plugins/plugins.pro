TEMPLATE = subdirs

SUBDIRS  = rational_function         \
           rational_data             \ 
			  rational_fitter_cgal      \
			  rational_fitter_quadprog  \
			  rational_fitter_quadproge \
			  rational_fitter_eigen     \
			  rational_fitter_leastsquare     \
			  rational_fitter_matlab

rational_fitter_cgal.depends      = rational_function rational_data
rational_fitter_quadprog.depends  = rational_function rational_data
rational_fitter_quadproge.depends = rational_function rational_data
rational_fitter_eigen.depends     = rational_function rational_data
rational_fitter_matlab.depends    = rational_function rational_data
