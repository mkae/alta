load(eigen)
load(ipopt)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, ipopt))


TARGET          = nonlinear_fitter_ipopt
TEMPLATE        = lib
CONFIG         *= plugin
DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore
