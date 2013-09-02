TARGET          = nonlinear_fitter_ceres
TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
						ceres   \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore


