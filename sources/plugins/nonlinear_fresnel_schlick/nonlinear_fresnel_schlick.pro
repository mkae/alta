TEMPLATE        = lib
CONFIG         *= plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = function.h
SOURCES         = function.cpp


LIBS           += -L../../build  \
						-lcore

