#pragma once

#include "function.h"
#include "data.h"
#include "args.h"
#include "common.h"

/*! \brief Fitting interface for generic fitting algorithms
 *  \ingroup core
 *
 *  \details
 */
class fitter
{
	public:

		// Virtual destructor
		virtual ~fitter() {}

		//! \brief static function to fit a data set d with the underling
		//! function class. Return the best fit (along with fitting
		//! information ?)
		virtual bool fit_data(const ptr<data> d, function* f, const arguments& args) = 0 ;

		//! \brief parse the command line arguments to setup some general 
		//! options before any fit. Those options should be resilient to 
		//! multiple call to the fit_data procedure
		virtual void set_parameters(const arguments& args) = 0 ;

} ;

