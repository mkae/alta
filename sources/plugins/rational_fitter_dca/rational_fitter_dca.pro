load(eigen)
load(matlab)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, matlab)) 

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
