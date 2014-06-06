#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/vertical_segment.h>

/*! \brief A non-linear fitter using the NLOpt solver
 *  \ingroup plugins
 *
 *  \details
 *  <h3>Third party requirements</h3>
 *  
 *  You will need three external libraries to compile this plugin:
 *  <ul>
 *		<li><a href="http://ab-initio.mit.edu/wiki/index.php/NLopt">NLopt</a> 
 *		library, version 2.3. On some linux distributions, you can install 
 *		libnlopt-dev using the package manager.</li>
 *  </ul>
 *
 *  You need to provide your own nlopt.prf file for qmake to generate the correct
 *  Makefile or Visual Studio solution. In particular this configuration file
 *  should provide:
 *
 *  <pre>
 *  INCLUDEPATH += [path-to-nlopt-include]
 *  LIBS += -L[path-to-nlopt-lib] -lnlopt -lm
 *  </pre>
 *
 *
 *  <h3>Plugin parameters</h3>
 *
 */
class nonlinear_fitter_nlopt: public fitter
{
	public: // methods

		nonlinear_fitter_nlopt() ;
		virtual ~nonlinear_fitter_nlopt() ;

		// Fitting a data object
		//
		virtual bool fit_data(const ptr<data> d, function* fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function


	protected: // data
} ;

