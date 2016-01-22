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

using namespace alta;

/*! \brief A non-linear fitter using the COIN-OR IpOpt solver
 *  \ingroup plugins
 *  \ingroup fitters
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
 *
 *  <h3>Plugin parameters</h3>
 *  <ul>
 *		<li><b>--ipopt-max-iter</b> <em>[int]</em> to control the number
 *		of iterations the non linear solver will take before returning a 
 *		solution. <em>Note that this solution might incorrect.</em></li>
 *		<li><b>--solver</b> <em>[string]</em> to control the type 
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
		virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function


	protected: // data
} ;

