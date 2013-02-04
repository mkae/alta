#pragma once 

#include <function.h>
#include <data.h>
#include <fitter.h>
#include <args.h>

class plugins_manager
{
	public: //functions

		// Create the object, parse the argument and load
		// all the plugins
		//
		plugins_manager(const arguments& args) ;

		// Get instances of the function, the data and the
		// fitter
		//
		function* get_function() const ;
		data*     get_data()     const ;
		fitter*   get_fitter()   const ;

	private: //data

		std::vector<function*> _functions ;
		std::vector<data*>     _datas ;
		std::vector<fitter*>   _fitters ;
} ;
