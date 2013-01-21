DEST_DIR = ../bin
CONFIG += debug
INCLUDEPATH += ../../ ../../libs/rational_1d

SOURCES += main.cpp ../../libs/rational_1d/rational_1d_fitter.cpp
QMAKE_CXXFLAGS += -std=c++11 -frounding-math

LIBS += -lCGAL -lboost_regex
