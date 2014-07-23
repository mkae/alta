load(matlab)
load(eigen)

requires(contains(CONFIG, eigen))
requires(contains(CONFIG, matlab)) 

TEMPLATE        = lib
CONFIG         *=  plugin 
DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp

LIBS           += -L../../build \
                  -lcore        

