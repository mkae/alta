load(openexr)
load(eigen)

requires(contains(CONFIG, openexr))
requires(contains(CONFIG, eigen))

TEMPLATE        = lib
CONFIG         *=  plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build   \
                  -lcore          
