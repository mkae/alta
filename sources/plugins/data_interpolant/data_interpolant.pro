TEMPLATE        = lib
CONFIG         *=  plugin
DESTDIR         = ../../build

load(flann)
load(eigen)

packagesExist(flann, eigen) {

	INCLUDEPATH    += ../..
	HEADERS         = data.h
	SOURCES         = data.cpp


	LIBS           += -L../../build \
							-lcore
}
