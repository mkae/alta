TARGET          = nonlinear_fitter_nlopt
TEMPLATE        = lib
CONFIG         *= plugin
DESTDIR         = ../../build

load(eigen)
load(nlopt)

packagesExist(eigen, nlopt) {
 
	INCLUDEPATH    += ../.. 

	HEADERS         = fitter.h
	SOURCES         = fitter.cpp

	LIBS           += -L../../build \
	                  -lcore
}

