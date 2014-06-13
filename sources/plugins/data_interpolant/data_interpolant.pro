TEMPLATE        = lib
CONFIG         *=  plugin	\
						 eigen	\
						 flann

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build \
						-lcore        \
						-lCGAL

unix {
	QMAKE_CXXFLAGS += -frounding-math 
}
