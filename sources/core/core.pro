TEMPLATE = lib
CONFIG  *= static  \
           console \
			  eigen

DESTDIR  = ../build

unix{
	QMAKE_CXXFLAGS += -std=c++0x -Wall -pedantic
}

HEADERS  = args.h               \
           common.h             \
		 	  data.h               \
	 		  fitter.h             \
	 		  function.h           \
 			  plugins_manager.h    \
			  vertical_segment.h   \
			  rational_function.h  \
			  params.h             \
           clustering.h

SOURCES  = common.cpp            \
           plugins_manager.cpp   \
           vertical_segment.cpp  \
		     rational_function.cpp \
		     params.cpp			   \
           function.cpp          \
#          clustering.cpp
