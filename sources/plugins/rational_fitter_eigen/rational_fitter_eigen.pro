load(eigen)

requires(contains(CONFIG, eigen))


TARGET          = rational_fitter_eigen
TEMPLATE        = lib
CONFIG         *= plugin

DESTDIR         = ../../build
 
INCLUDEPATH    += ../rational_function \
                  ../rational_data     \
                  ../.. 

HEADERS         = rational_fitter.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build        \
                  -lcore
