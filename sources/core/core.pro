TEMPLATE = lib
CONFIG  *= static  \
           qt      \
           console

DESTDIR  = ../build

unix{
	QMAKE_CXXFLAGS += -std=c++0x -m64
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

SOURCES  = plugins_manager.cpp   \
           vertical_segment.cpp  \
		   rational_function.cpp \
		   params.cpp			 \
           function.cpp          \
#          clustering.cpp
    rational_function.inl
