TARGET          = nonlinear_fitter_nlopt
TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
						nlopt   \
						eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore


