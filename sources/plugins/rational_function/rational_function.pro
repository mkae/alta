TEMPLATE        = lib
CONFIG         += plugin 
 
INCLUDEPATH    += ../..
HEADERS         = rational_function.h
SOURCES         = rational_function.cpp

LIBS           += -lboost_regex

QMAKE_CXXFLAGS += -std=c++11 -frounding-math 

