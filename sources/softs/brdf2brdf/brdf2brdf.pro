CONFIG         += console \
                  eigen
CONFIG         -= app_bundle


DESTDIR         = ../../build
INCLUDEPATH    += ../../ 

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

unix{
   PRE_TARGETDEPS += ../../build/libcore.a
	LIBS += -ldl
}
