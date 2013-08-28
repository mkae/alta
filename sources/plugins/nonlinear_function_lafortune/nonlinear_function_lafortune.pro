TEMPLATE        = lib
CONFIG         *=  qt     \
                   plugin \
						 eigen

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = function.h
SOURCES         = function.cpp


LIBS           += -L../../build  \
						-lcore

