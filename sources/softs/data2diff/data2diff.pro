CONFIG         += qt	eigen

INCLUDEPATH    += ../../
DESTDIR         = ../../build

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

unix{
   PRE_TARGETDEPS += ../../build/libcore.a
	LIBS += -ldl
}
