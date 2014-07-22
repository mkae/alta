TARGET          = nonlinear_fitter_ipopt
TEMPLATE        = lib
CONFIG         *= plugin
DESTDIR         = ../../build

load(eigen)
load(ipopt)

packagesExist(eigen, ipopt) {

 
	INCLUDEPATH    += ../.. 

	HEADERS         = fitter.h
	SOURCES         = fitter.cpp

	LIBS           += -L../../build \
	                  -lcore
}

