load(eigen)
load(quadprog)
load(parsec)
load(coreblas)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, quadprog))
requires(contains(CONFIG, parsec))
requires(contains(CONFIG, coreblas))


TARGET    = rational_fitter_parsec_one
TEMPLATE  = lib
CONFIG   *= plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../rational_function \
	               ../rational_data     \
		            ../..

HEADERS         = rational_fitter.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build        \
                  -lcore

QMAKE_CXXFLAGS += -g3 -O0
