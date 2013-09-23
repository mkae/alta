DEST_DIR = ../bin
CONFIG += plugin	\
#         debug
QT += 
INCLUDEPATH += ../../ ../../plugins/rational_1d /home/belcour/Sources/Eigen/include/eigen3

SOURCES += main.cpp															\
           ../../plugins/rational_1d/rational_1d_fitter.cpp       \
           ../../plugins/rational_1d/rational_1d_fitter_cgal.cpp  \
           ../../plugins/rational_1d/rational_1d_fitter_eigen.cpp

QMAKE_CXXFLAGS += -std=c++11 -frounding-math 

LIBS += -lCGAL -lboost_regex
LIBS += -ldl
