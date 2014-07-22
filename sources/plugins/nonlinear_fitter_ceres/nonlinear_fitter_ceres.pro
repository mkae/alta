TARGET          = nonlinear_fitter_ceres
TEMPLATE        = lib
CONFIG         *= plugin

load(ceres)
load(eigen)

packagesExist(ceres, eigen) {
	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../.. 

	HEADERS         = fitter.h
	SOURCES         = fitter.cpp

	LIBS           += -L../../build \
		               -lcore
}

