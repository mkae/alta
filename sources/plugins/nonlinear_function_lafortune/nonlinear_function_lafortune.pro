TEMPLATE        = lib
CONFIG         *=  qt     \
                   plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = function.h
SOURCES         = function.cpp


LIBS           += -L../../build  \
						-lcore

