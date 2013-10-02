TARGET          = rational_fitter_quadprog
TEMPLATE        = lib
CONFIG         *= plugin     \
						eigen      \
						quadprog 

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build        \
						-lcore	

#QMAKE_CXXFLAGS += -fPIC

