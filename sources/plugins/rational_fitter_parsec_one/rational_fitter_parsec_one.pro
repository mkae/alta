TARGET    = rational_fitter_parsec_one
TEMPLATE  = lib
CONFIG   *= plugin

load(eigen)
load(quadprog)
load(parsec)
load(coreblas)

packagesExist(eigen, quadprog, parsec, coreblas) {

	DESTDIR         = ../../build

	INCLUDEPATH    += ../rational_function \
		               ../rational_data     \
			            ../..

	HEADERS         = rational_fitter.h
	SOURCES         = rational_fitter.cpp

	LIBS           += -L../../build        \
	                  -lcore
	
	QMAKE_CXXFLAGS += -g3 -O0
}
