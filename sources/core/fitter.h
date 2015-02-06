/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

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
		virtual bool fit_data(const ptr<data>& d, ptr<function>& f, const arguments& args) = 0 ;

		//! \brief parse the command line arguments to setup some general
		//! options before any fit. Those options should be resilient to
		//! multiple call to the fit_data procedure
		virtual void set_parameters(const arguments& args) = 0 ;

} ;
