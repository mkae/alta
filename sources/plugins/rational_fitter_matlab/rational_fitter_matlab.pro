TARGET          = rational_fitter_matlab
TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
						eigen   \
						matlab  

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build           \
                  -lcore

