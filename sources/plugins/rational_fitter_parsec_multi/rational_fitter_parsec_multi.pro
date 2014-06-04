TARGET   = rational_fitter_parsec_multi
TEMPLATE = lib
CONFIG	 *= plugin      \
	eigen           \
	quadprog        \
	parsec          \
	coreblas

DESTDIR         = ../../build

INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../..

HEADERS         = rational_fitter.h   \
                  quadratic_program.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build       \
		-lcore

#QMAKE_CXXFLAGS += -fPIC
QMAKE_CXXFLAGS += -fPIC -g

