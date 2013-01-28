DEST_DIR = ../bin
CONFIG += debug plugin
QT += 
INCLUDEPATH += ../../ ../../libs/rational_1d /home/belcour/Sources/Eigen/include/eigen3

SOURCES += main.cpp

QMAKE_CXXFLAGS += -std=c++11 -frounding-math -fPIC
QMAKE_LFLAGS   +=  -Wl,-rpath="/home/belcour/Projects/alta/sources/plugins/build"

LIBS += -lCGAL -lboost_regex
