TARGET          = rational_fitter_dca
TEMPLATE        = lib

load(eigen)
load(matlab)

packagesExist(eigen, matlab) {
	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../rational_function \
	                  ../rational_data     \
	                  ../.. 

	HEADERS         = rational_fitter.h
	SOURCES         = rational_fitter.cpp

	LIBS           += -L../../build \
	                  -lcore
}
