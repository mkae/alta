TEMPLATE = lib
CONFIG  *= static \
           qt

DESTDIR  = ../build

HEADERS  = args.h               \
           common.h             \
		 	  data.h               \
	 		  fitter.h             \
	 		  function.h           \
 			  plugins_manager.h    \
			  vertical_segment.h   \
			  rational_function.h  \
			  params.h

SOURCES  = plugins_manager.cpp  \
           vertical_segment.cpp \
			  rational_function.cpp
