TEMPLATE        = lib
CONFIG         *=  plugin 

load(matlab)
load(eigen)

packagesExist(matlab, eigen) {
	DESTDIR         = ../../build

	INCLUDEPATH    += ../..
	HEADERS         = data.h
	SOURCES         = data.cpp


	LIBS           += -L../../build \
							-lcore        
}
