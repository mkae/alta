DEST_DIR = ../bin
CONFIG += debug
INCLUDEPATH += ../../ ../../libs/rational_1d /home/belcour/Sources/Eigen/include/eigen3

SOURCES += main.cpp															\
           ../../libs/rational_1d/rational_1d_fitter.cpp       \
           ../../libs/rational_1d/rational_1d_fitter_cgal.cpp  \
           ../../libs/rational_1d/rational_1d_fitter_eigen.cpp

QMAKE_CXXFLAGS += -std=c++11 -frounding-math 

LIBS += -lCGAL -lboost_regex
