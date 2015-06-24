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

/*! \brief A fitter for non-linear BRDF models that uses Eigen's
 *  Levenberg-Marquardt solver.
 *  \ingroup plugins
 *  \ingroup fitters
 *
 *  \details
 *  You will need the Eigen unsupported library (available with mercurial)
 *
 *  #### Plugin parameters
 *
 *	 + `--fit-compound` to control how the fitting procedure is done. If this
 *	 flag is set, any compound function will be decomposed during the fit. The
 *	 fitting will be done incrementally. *in progress*
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

