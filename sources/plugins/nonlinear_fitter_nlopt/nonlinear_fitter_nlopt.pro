load(eigen)
load(nlopt)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, nlopt))


TARGET          = nonlinear_fitter_nlopt
TEMPLATE        = lib
CONFIG         *= plugin
DESTDIR         = ../../build

INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
                  -lcore