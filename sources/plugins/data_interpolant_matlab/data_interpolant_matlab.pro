TEMPLATE        = lib
CONFIG         *=  plugin 
DESTDIR         = ../../build

load(matlab)
load(eigen)

!contains(CONFIG, eigen)|!contains(CONFIG, matlab) {
	unset(TARGET)
} else {

	INCLUDEPATH    += ../..
	HEADERS         = data.h
	SOURCES         = data.cpp


	LIBS           += -L../../build \
							-lcore        
}
