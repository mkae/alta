TEMPLATE        = lib
CONFIG         *=  plugin

load(flann)
load(eigen)

packagesExist(flann, eigen) {
	DESTDIR         = ../../build

	INCLUDEPATH    += ../..
	HEADERS         = data.h
	SOURCES         = data.cpp


	LIBS           += -L../../build \
							-lcore
}
