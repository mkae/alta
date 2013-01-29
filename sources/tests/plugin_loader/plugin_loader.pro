DESTDIR = ../../build
CONFIG += debug plugin
QT += 
INCLUDEPATH += ../../ ../../libs/rational_1d \
#					/home/belcour/Sources/Eigen/include/eigen3

SOURCES += main.cpp

QMAKE_CXXFLAGS += -frounding-math -fPIC
QMAKE_LFLAGS   +=  -Wl,-rpath='\$\$ORIGIN:.:./build:./plugins'

LIBS += -lCGAL
