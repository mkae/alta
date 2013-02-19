TEMPLATE        = lib
CONFIG         *=  qt     \
                   plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build \
						-lcore


