load(eigen)

requires(contains(CONFIG, eigen))

TARGET          = nonlinear_fitter_eigen
TEMPLATE        = lib
CONFIG         *= plugin

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore


