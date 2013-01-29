TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin
						
DESTDIR         = ../../build
 
INCLUDEPATH    += ../..
HEADERS         = rational_data.h
SOURCES         = rational_data.cpp

QMAKE_CXXFLAGS += -frounding-math -fPIC
