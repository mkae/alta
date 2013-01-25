TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin
						
DESTDIR         = ../build
 
INCLUDEPATH    += ../..
HEADERS         = rational_data.h
SOURCES         = rational_data.cpp

LIBS           += -lboost_regex

QMAKE_CXXFLAGS += -std=c++11 -frounding-math -fPIC
