TEMPLATE = lib
CONFIG  *= static \
           qt

DESTDIR  = ../build

HEADERS  = args.h              \
           common.h            \
		 	  data.h              \
	 		  fitter.h            \
	 		  function.h          \
 			  plugin_manager.h

SOURCES  = plugin_manager.cpp
