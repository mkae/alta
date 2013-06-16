CONFIG         += qt \
                  console

DESTDIR         = ../../build
INCLUDEPATH    += ../../ ../../libs/rational_1d \

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

#QMAKE_CXXFLAGS += -frounding-math -fPIC
#QMAKE_LFLAGS   +=  -Wl,-rpath='\$\$ORIGIN:.:./build:./plugins'
