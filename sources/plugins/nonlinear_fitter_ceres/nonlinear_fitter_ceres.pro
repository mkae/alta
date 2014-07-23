load(ceres)
load(eigen)

requires(contains(CONFIG, ceres)) 
requires(contains(CONFIG, eigen))


TARGET          = nonlinear_fitter_ceres
TEMPLATE        = lib
CONFIG         *= plugin
DESTDIR         = ../../build

INCLUDEPATH    += ../.. 

HEADERS         = fitter.h
SOURCES         = fitter.cpp

LIBS           += -L../../build \
	               -lcore
