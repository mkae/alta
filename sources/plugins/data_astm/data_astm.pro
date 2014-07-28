TEMPLATE        = lib
CONFIG         *= plugin
DESTDIR         = ../../build

load(eigen)

packagesExist(ceres, eigen) {

	INCLUDEPATH    += ../..
	HEADERS         = data.h
	SOURCES         = data.cpp


	LIBS           += -L../../build \
							-lcore
}
