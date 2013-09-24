CONFIG         += qt      \
                  console \
						eigen

DESTDIR         = ../../build
INCLUDEPATH    += ../../ 

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

unix{
   PRE_TARGETDEPS += ../../build/libcore.a
	LIBS += -ldl
}
