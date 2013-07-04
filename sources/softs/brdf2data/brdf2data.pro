CONFIG         += qt      \
                  console

INCLUDEPATH    += ../../
DESTDIR         = ../../build

SOURCES        += main.cpp
LIBS           += -L../../build -lcore
