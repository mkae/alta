TEMPLATE        = lib
CONFIG         += plugin  \
						eigen   \
						cgal    
CONFIG         -= qt

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = rational_fitter_cgal.h
SOURCES         = rational_fitter_cgal.cpp

LIBS           += -L../../build        \
						-lcore

unix {
	QMAKE_CXXFLAGS += -frounding-math 
}
