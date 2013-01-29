TEMPLATE        = lib
CONFIG         *= qt     \
                  plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = rational_function.h
SOURCES         = rational_function.cpp

#LIBS           += -lboost_regex

QMAKE_CXXFLAGS += -frounding-math -fPIC -rdynamic -g

