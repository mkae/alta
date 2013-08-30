TEMPLATE        = lib
CONFIG         *=  qt     \
                   plugin \
						 eigen

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build \
						-lcore


