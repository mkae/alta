TEMPLATE        = lib
CONFIG         *= qt     \
                  plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = rational_function.h
SOURCES         = rational_function.cpp

QMAKE_CXXFLAGS += -frounding-math \
                  -fPIC           \
						-g

