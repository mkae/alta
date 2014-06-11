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

/*! \brief A fitter for non-linear BRDF models that uses Eigen's
 *  Levenberg-Marquardt solver.
 *  \ingroup plugins
 *
 * You will need the Eigen unsupported library (available with mercurial)
 *  <ul>
 *		<li>Eigen latest version.</li>
 *  </ul>
 *
 *  You need to provide your own eigen.prf file for qmake to generate the correct
 *  Makefile or Visual Studio solution. In particular this configuration file
 *  should provide:
 *
 *  <pre>
 *  INCLUDEPATH += [path-to-eigen-include]
 *  </pre>
 *
 *
 *  <h3>Plugin parameters</h3>
 *  <ul>
 *		<li><b>--fit-compound</b> to control how the fitting procedure is
 *		done. If this flag is set, any compound function will be decomposed
 *      during the fit. The fitting will be done incrementally. <b>in progress</b></li>
 *  </ul>
 */
class nonlinear_fitter_eigen: public fitter
{
	public: // methods

		nonlinear_fitter_eigen() ;
		virtual ~nonlinear_fitter_eigen() ;

		// Fitting a data object
		//
		virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function


	protected: // data
} ;

