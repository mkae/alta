TARGET          = rational_fitter_leastsquare

TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build           \
                  -lcore


#QMAKE_CXXFLAGS += -fPIC


