/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

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
 *  \ingroup fitters
 *
 *  \details
 *  #### Third party requirements
 *  
 *  You will need three external libraries to compile this plugin:
 *   
 *   + The [NlOpt library, version 2.3][nlopt]. On some linux distributions, you
 *     can install libnlopt-dev using the package manager.
 *   + If NlOpt is not installed on your system, it will be downloaded in the
 *     `external` directory and compiled on Unix systems (GNU/Linux, OSX, ...).
 *     However, this compilation requires a Fortran compiler.
 *
 *
 *  #### Plugin parameters
 *
 *   + `--nlopt-optimizer [string]` permits to select one of [NlOpt's various
 *     optimizer][optimizers]. By default it selects the Local Sequential
 *     Quadratic Programming algorithm (NLOPT_LD_SLSQP).
 *
 *   + `--nlop-max-num-iterations [int]` permits to change the number of iterations
 *     of NlOpt before stopping and returning a result. By default, the value is
 *     `10`.
 *   + '--nlop-relative-function-tolerance [float]'. Default value is 1e-4.
 *   + '--nlop-abs-function-tolerance [float]'. Default valie is 1e-6.
 *
 *  [nlopt]: http://ab-initio.mit.edu/wiki/index.php/NLopt
 *  [optimizers]: http://ab-initio.mit.edu/wiki/index.php/NLopt
 */
class nonlinear_fitter_nlopt: public fitter
{
	public: // methods

		nonlinear_fitter_nlopt() ;
		virtual ~nonlinear_fitter_nlopt() ;

		// Fitting a data object
		//
		virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function


	protected: // data
} ;

