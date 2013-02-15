TARGET          = rational_fitter_quadproge
TEMPLATE        = lib
CONFIG         *= qt         \
                  plugin     \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h   \
                  eiquadprog.hpp
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build        \
                  -lcore

