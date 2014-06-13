TARGET          = nonlinear_fitter_nlopt
TEMPLATE        = lib
CONFIG         *= plugin

load(eigen)
load(nlopt)

packagesExist(eigen, nlopt) {
	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../.. 

	HEADERS         = fitter.h
	SOURCES         = fitter.cpp

	LIBS           += -L../../build \
	                  -lcore
}

