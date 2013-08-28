CONFIG         += qt      \
                  console \
						eigen

DESTDIR         = ../../build
INCLUDEPATH    += ../../ ../../libs/rational_1d \

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

unix{
   PRE_TARGETDEPS += ../../build/libcore.a

	CONFIG += debug
}
