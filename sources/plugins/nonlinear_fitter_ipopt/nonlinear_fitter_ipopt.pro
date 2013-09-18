TARGET          = nonlinear_fitter_ipopt
TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
						ipopt   \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore


