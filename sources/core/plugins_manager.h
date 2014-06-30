#pragma once

#include <map>
#include <string>

#include "args.h"
#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"
#include "clustering.h"
#include "ptr.h"

/*! \class plugins_manager
 *  \brief This class permits to load plugin from shared library files.
 *  \ingroup core
 *
 *  \details
 *  This class handles the loading of plugins and insure that they can
 *  talk to each others through coordinates transforms.
 */
class plugins_manager
{
	public: //functions

		//! \brief get an instance of the function that is defined in the plugin with
		//! filename n. Return null if no one exist.
		//!
		//! \details
		//! This function try to load the shared object file specified in the
		//! <code>--func filename</code>.
		static function* get_function(const arguments& args) ;

		//! \brief load a function from the ALTA input file.
		static function* get_function(const std::string& filename);

		//! \brief get an instance of the data that is defined in the plugin with
		//! filename n. Return null if no one exist.
		static ptr<data> get_data(const std::string& n) ;

		//! \brief get an instance of the fitter that is defined in the plugin with
		//! filename n. Return null if no one exist.
		static ptr<fitter> get_fitter(const std::string& n) ;


		//! \brief check if a data object and a function object are compatibles.
		//! this has to be done before fitting to ensure that the
		//! parametrizations spaces are the same.
		//! \todo specify an output parametrization for the function ?
		static void check_compatibility(ptr<data>& d, const ptr<function>& f,
				const arguments& args) ;


		//! \brief Provide a measure of how much memory there is on the system.
		//! \details It permits to know is one can allocate more memory for a fitting
		//! procedure for example.
		static size_t get_system_memory() ;

	private: //data

		// Object provider prototypes
		typedef function* (*FunctionPrototype)();
		typedef fitter*   (*FitterPrototype)();
		typedef data*     (*DataPrototype)();
} ;
