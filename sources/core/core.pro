TEMPLATE = lib
CONFIG  *= static \
           qt

DESTDIR  = ../build

HEADERS  = args.h              \
           common.h            \
		 	  data.h              \
	 		  fitter.h            \
	 		  function.h          \
 			  plugins_manager.h

SOURCES  = plugins_manager.cpp
