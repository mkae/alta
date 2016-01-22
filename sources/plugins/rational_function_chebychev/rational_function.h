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
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

using namespace alta;

class rational_function_chebychev_1d : public rational_function_1d
{
public: // methods

    rational_function_chebychev_1d() ;
    rational_function_chebychev_1d(int nX, int np, int nq) ;
    virtual ~rational_function_chebychev_1d() {}

    // Get the p_i and q_j function
    virtual double p(const vec& x, int i) const ;
    virtual double q(const vec& x, int j) const ;

protected:  // methods

} ;

/*! \ingroup functions
 *  \class rational_function_chebychev
 *  \brief Rational function using 
 *  [Chebychev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials).
 *
 *  \details 
 *  This rational function uses Chebychev polynomials as the basis function for the
 *  numerator and the denominator of the rational function: 
 *  <center>
 *  \f$ f(x) = \sum a_i p_i(x) / b_i q_{i}(x) \f$, with \f$ p_i = q_i\f$.
 *  </center>
 *
 *  Chebychev polynomials can be defined using trigonometric functions:
 *  \f$ p_i(x) = \cos\left( i \; \mbox{acos}(x) \right)\f$. We use this formulation.
 *
 *  For more details see the addendum on Rational BRDF available at https://hal.inria.fr/hal-00913516
 *
 *  \author Laurent Belcour \<laurent.belcour@umontreal.ca\>
 */
class rational_function_chebychev : public rational_function
{
	public: // methods

		rational_function_chebychev() ;
		virtual ~rational_function_chebychev() ;
		
		//! Get the 1D function associated with color channel i. If no one exist, 
		//! this function allocates a new element. If i > nY, it returns NULL.
		virtual rational_function_1d* get(int i) ;

		//! Update the y-1D function for the ith dimension.
		//! \note It will test if the 1D function provided is of the dynamic type
		//! \name rational_function_chebychev_1d
		virtual void update(int i, rational_function_1d* r)
		{
			if(dynamic_cast<rational_function_chebychev_1d*>(r) != NULL)
			{
				rational_function::update(i, r);
			}
			else
			{
#ifdef DEBUG
				std::cerr << "<<ERROR>> the function provided is not of type \"rational_function_chebychev\"" << std::endl;
#endif
			}
		}
} ;

