#pragma once 

#include <map>
#include <string>

#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"

class plugins_manager
{
	public: //functions

		// Create the object, parse the argument and load
		// all the plugins
		//
		plugins_manager(const arguments& args) ;

		// Get instances of the function, the data and the
		// fitter. Select the first in the map,
		//
		function* get_function() const ;
		data*     get_data()     const ;
		fitter*   get_fitter()   const ;
		
		// Get instances of the function, the data and the
		// fitter, select one based on the name. Return null
		// if no one exist.
		//
		function* get_function(const std::string& n) const ;
		data*     get_data(const std::string& n)     const ;
		fitter*   get_fitter(const std::string& n)   const ;

	private: //data

		std::map<std::string, function*> _functions ;
		std::map<std::string, data*>     _datas ;
		std::map<std::string, fitter*>   _fitters ;
} ;
