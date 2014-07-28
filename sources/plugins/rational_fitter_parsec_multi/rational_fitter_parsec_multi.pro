load(eigen)
load(quadprog)
load(parsec)
load(coreblas)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, quadprog))
requires(contains(CONFIG, parsec))
requires(contains(CONFIG, coreblas))


TARGET    = rational_fitter_parsec_multi
TEMPLATE  = lib
CONFIG	 *= plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../..

HEADERS         = rational_fitter.h   \
                  quadratic_program.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build       \
                  -lcore

!debug:QMAKE_CXXFLAGS += -fPIC
debug:QMAKE_CXXFLAGS  += -fPIC -g