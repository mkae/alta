CONFIG         += qt console eigen

INCLUDEPATH    += ../../ ../../plugins
DESTDIR         = ../../build

SOURCES        += main.cpp \
						../../plugins/nonlinear_fresnel_schlick/function.cpp
LIBS           += -L../../build -lcore

unix{
   PRE_TARGETDEPS += ../../build/libcore.a
}
