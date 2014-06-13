TEMPLATE        = lib
CONFIG         += plugin
CONFIG         -= qt

load(eigen)
load(cgal)

packagesExist(eigen, cgal) {
	DESTDIR         = ../../build
 
	INCLUDEPATH    += ../.. 

	HEADERS         = rational_fitter_cgal.h
	SOURCES         = rational_fitter_cgal.cpp

	LIBS           += -L../../build        \
							-lcore

	unix {
		QMAKE_CXXFLAGS += -frounding-math 
	}
}
