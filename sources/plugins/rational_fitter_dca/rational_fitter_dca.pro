
load(eigen)
load(matlab)

!contains(CONFIG, eigen)|!contains(CONFIG, matlab) {
	unset(TARGET)
} else {
	TARGET          = rational_fitter_dca
	TEMPLATE        = lib
	CONFIG         *= plugin
	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../rational_function \
	                  ../rational_data     \
	                  ../.. 

	HEADERS         = rational_fitter.h
	SOURCES         = rational_fitter.cpp

	LIBS           += -L../../build \
	                  -lcore
}
