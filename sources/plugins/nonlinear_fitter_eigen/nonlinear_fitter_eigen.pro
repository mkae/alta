TARGET          = nonlinear_fitter_eigen
TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore

