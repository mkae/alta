TARGET          = quadprog++
TEMPLATE        = lib
CONFIG         *= static  \
				  qt      \
                  plugin  \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = QuadProg++.hh \
                  Array.hh
SOURCES         = QuadProg++.cc \
                  Array.cc

