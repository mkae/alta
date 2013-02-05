CONFIG         += qt

DESTDIR         = ../../build
INCLUDEPATH    += ../../ ../../libs/rational_1d \

SOURCES        += main.cpp

QMAKE_CXXFLAGS += -frounding-math -fPIC
#QMAKE_LFLAGS   +=  -Wl,-rpath='\$\$ORIGIN:.:./build:./plugins'
