CONFIG         += qt     \
						matlab  

INCLUDEPATH    += ../../
DESTDIR         = ../../build

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

unix{
	LIBS += -ldl
}
