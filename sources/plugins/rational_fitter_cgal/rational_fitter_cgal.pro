TEMPLATE        = lib
CONFIG         *= shared  \
                  qt      \
                  plugin  \
						eigen   \
						cgal    

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter_cgal.h
SOURCES         = rational_fitter_cgal.cpp

LIBS           += -L../../build           \
						-lrational_function	\
						-lrational_data

QMAKE_CXXFLAGS += -frounding-math -fPIC

