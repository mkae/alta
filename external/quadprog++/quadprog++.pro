TARGET          = quadprog++
TEMPLATE        = lib
CONFIG         *= qt      \
                  plugin  \
                  eigen

DESTDIR         = ../../build
 
INCLUDEPATH    += ../.. 

HEADERS         = QuadProg++.hh \
                  Array.hh
SOURCES         = QuadProg++.cc \
                  Array.cc

