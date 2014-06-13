TARGET          = nonlinear_fitter_ipopt
TEMPLATE        = lib
CONFIG         *= plugin

load(eigen)
load(ipopt)

packagesExist(eigen, ipopt) {

	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../.. 

	HEADERS         = fitter.h
	SOURCES         = fitter.cpp

	LIBS           += -L../../build \
	                  -lcore
}

