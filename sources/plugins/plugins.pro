TEMPLATE = subdirs
SUBDIRS  = rational_function    \
           rational_data        \ 
			  rational_fitter_cgal 

rational_fitter_cgal.depends = rational_function rational_data
