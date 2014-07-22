
load(cgal)
load(eigen)


!contains(CONFIG, eigen)|!contains(CONFIG, cgal) {
	unset(TARGET)
} else {
	TARGET          = rational_fitter_cgal
	TEMPLATE        = lib
	CONFIG         += plugin
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
}
