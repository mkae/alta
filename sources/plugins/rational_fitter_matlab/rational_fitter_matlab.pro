load(eigen)
load(matlab)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, matlab)) 


TARGET          = rational_fitter_matlab
TEMPLATE        = lib
CONFIG         *= plugin
CONFIG         -= qt

DESTDIR         = ../../build

INCLUDEPATH    += ../.. 

HEADERS         = rational_fitter.h
SOURCES         = rational_fitter.cpp

LIBS           += -L../../build        \
	               -lcore