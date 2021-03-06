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
#include <engine.h>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/rational_function.h>
#include <core/vertical_segment.h>

using namespace alta;

/*! \brief A rational function optimizer following the DCA algorithm.
 * \ingroup plugins
 * \ingroup fitters
 * \todo Implement Papamarkos fitter?
 * \todo I should be able to test when load a BRDF text file to ensure the
 * loaded object is correct.
 */
class rational_fitter_dca : public fitter
{

	public: // methods
	
		rational_fitter_dca() ;
		virtual ~rational_fitter_dca() ;
			
		// Fitting a data object
		//
        virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
        virtual bool fit_data(const ptr<data>& d, ptr<rational_function>& fit, const arguments& args) ;

        //! \brief Bootstrap the DCA algorithm with an already fitted function. It will
		//! load the the rational function object from a text file defined in the argument
		//! --bootstrap %filename%.
		void bootstrap(const ptr<data> d, int& np, int& nq, 
		               const ptr<rational_function>& fit, double& delta,
		               const arguments& args) ;

	protected: // data

        // Matlab Engine
		Engine *ep;
} ;

