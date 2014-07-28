load(eigen)
load(quadprog)
load(openmp)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, quadprog))


TARGET          = rational_fitter_parallel
TEMPLATE        = lib
CONFIG         *= plugin

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h   \
                  quadratic_program.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build       \
                  -lcore	