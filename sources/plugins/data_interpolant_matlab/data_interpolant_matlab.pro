TEMPLATE        = lib
CONFIG         *=  plugin \
						 eigen  \
                   matlab  

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build \
						-lcore        

