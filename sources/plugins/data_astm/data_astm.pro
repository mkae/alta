TEMPLATE        = lib
CONFIG         *= plugin


load(eigen)

packagesExist(ceres, eigen) {
	DESTDIR         = ../../build

	INCLUDEPATH    += ../..
	HEADERS         = data.h
	SOURCES         = data.cpp


	LIBS           += -L../../build \
							-lcore
}
