TEMPLATE        = lib
CONFIG         *=  plugin  \
						 eigen   \
						 openexr

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build   \
						-lcore          
