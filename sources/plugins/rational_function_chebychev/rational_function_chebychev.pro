TEMPLATE        = lib
CONFIG         *=  qt     \
                   plugin

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = rational_function.h
SOURCES         = rational_function.cpp


LIBS           += -L../../build        \
						-lcore

#QMAKE_CXXFLAGS += -frounding-math \
#                  -fPIC           \
#						-g

