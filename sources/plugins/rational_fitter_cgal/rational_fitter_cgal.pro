TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
						debug

DESTDIR         = ../build
 
INCLUDEPATH    += ../rational_function ../rational_data ../.. /home/belcour/Sources/Eigen/include/eigen3
HEADERS         = rational_1d_fitter_cgal.h
SOURCES         = rational_1d_fitter_cgal.cpp

LIBS           += -lCGAL -lboost_regex                 \
                  -L../build                           \
						-lrational_function -lrational_data

QMAKE_CXXFLAGS += -std=c++11 -frounding-math -fPIC -g

