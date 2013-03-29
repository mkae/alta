TARGET          = rational_fitter_parallel
TEMPLATE        = lib
CONFIG         *= qt         \
                  plugin     \
						eigen      \
						quadprog   \
						openmp

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h \
    quadratic_program.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build        \
						-lcore	

#QMAKE_CXXFLAGS += -fPIC

