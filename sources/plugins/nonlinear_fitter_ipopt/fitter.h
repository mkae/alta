#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>

/*! \brief A non-linear fitter using the COIN-OR IpOpt solver
 *  \ingroup plugins
 *
 *  \details
 *  <h3>Third party requirements</h3>
 *  
 *  You will need three external libraries to compile this plugin:
 *  <ul>
 *		<li><a href="http://www.coin-or.org/Ipopt/">IpOpt</a> library, 
 *		version 3.11.3 and all its dependancies</li>
 *  </ul>
 *
 *  You need to provide your own ipopt.prf file for qmake to generate the correct
 *  Makefile or Visual Studio solution. In particular this configuration file
 *  should provide:
 *
 *  <pre>
 *  INCLUDEPATH += [path-to-ipopt-include]
 *  LIBS += -L[path-to-ipopt-lib] -lipopt -lipoptamplinterface 
 *  -lcoinmumps -lcoinmetis -lcoinlapack -lcoinblas -lcoinasl 
 *  </pre>
 *
 *
 *  <h3>Plugin parameters</h3>
 *  <ul>
 *		<li><b>--ipopt-max-iter</b> <em>[int]</em> to control the number
 *		of iterations the non linear solver will take before returning a 
 *		solution. <em>Note that this solution might incorrect.</em></li>
 *		<li><b>--ceres-solverr</b> <em>[string]</em> to control the type 
 *		linear solver that will be used during matrix operations. See
 *		<a href="http://www.coin-or.org/Ipopt/documentation/node50.html">
 *		here</a> for the list of possible choices.</li>
 *  </ul>
 *
 */
class nonlinear_fitter_ipopt: public fitter
{
	public: // methods

		nonlinear_fitter_ipopt() ;
		virtual ~nonlinear_fitter_ipopt() ;

		// Fitting a data object
		//
		virtual bool fit_data(const ptr<data> d, function* fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function


	protected: // data
} ;

