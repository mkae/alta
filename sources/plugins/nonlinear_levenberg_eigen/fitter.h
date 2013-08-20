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
 */
class nonlinear_fitter_eigen: public fitter
{
	public: // methods

		nonlinear_fitter_eigen() ;
		virtual ~nonlinear_fitter_eigen() ;

		// Fitting a data object
		//
		virtual bool fit_data(const data* d, function* fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function


	protected: // data
} ;

