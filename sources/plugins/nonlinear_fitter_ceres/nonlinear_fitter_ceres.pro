
load(ceres)
load(eigen)

!contains(CONFIG, ceres)|!contains(CONFIG, eigen) {
	unset(TARGET)
} else {
	TARGET          = nonlinear_fitter_ceres
	TEMPLATE        = lib
	CONFIG         *= plugin
	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../.. 

	HEADERS         = fitter.h
	SOURCES         = fitter.cpp

	LIBS           += -L../../build \
		               -lcore
}

