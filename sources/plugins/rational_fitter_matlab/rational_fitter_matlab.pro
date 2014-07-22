
load(eigen)
load(matlab)

!contains(CONFIG, eigen)|!contains(CONFIG, matlab) {
	unset(TARGET)
} else {
	TARGET          = rational_fitter_matlab
	TEMPLATE        = lib
	CONFIG         *= plugin  
	DESTDIR         = ../../build

	INCLUDEPATH    += ../.. 

	HEADERS         = rational_fitter.h
	SOURCES         = rational_fitter.cpp

	LIBS           += -L../../build        \
		               -lcore
}
