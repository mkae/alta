TEMPLATE        = lib
CONFIG         += plugin 
 
INCLUDEPATH    += ../rational_function ../rational_data ../.. /home/belcour/Sources/Eigen/include/eigen3
HEADERS         = rational_1d_fitter_cgal.h
SOURCES         = rational_1d_fitter_cgal.cpp

LIBS           += -lCGAL -lboost_regex                       \
                  -L../rational_function -L../rational_data  \
						-lrational_function -lrational_data

QMAKE_CXXFLAGS += -std=c++11 -frounding-math 

