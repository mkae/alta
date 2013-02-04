TEMPLATE        = lib
CONFIG         *= static  \
                  qt      \
						
DESTDIR         = ../../build
 
INCLUDEPATH    += ../..
HEADERS         = rational_data.h
SOURCES         = rational_data.cpp

QMAKE_CXXFLAGS += -frounding-math -fPIC
