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

class rational_function_legendre_1d : public rational_function_1d
{
	public: // methods

		rational_function_legendre_1d() ;
		rational_function_legendre_1d(int nX, int np, int nq) ;
		virtual ~rational_function_legendre_1d() {}

		// Get the p_i and q_j function
		virtual double p(const vec& x, int i) const ;
		virtual double q(const vec& x, int j) const ;

	protected:  // methods

		// Legendre polynomial evaluation
		double legendre(double x, int i) const;
} ;

/*! \ingroup functions
 *  \brief Rational function using Legendre polynomials.
 *  Plugin name: *rational_function_legendre*
 *
 *  \details 
 *
 *  \author Laurent Belcour <laurent.belcour@gmail.com>
 *
 */
class rational_function_legendre : public rational_function
{
	public: // methods

		rational_function_legendre() ;
		virtual ~rational_function_legendre() ;

		//! Get the 1D function associated with color channel i. If no one exist, 
		//! this function allocates a new element. If i > nY, it returns NULL.
		virtual rational_function_1d* get(int i)
		{
			// Check for consistency in the index of color channel
			if(i < _nY)
			{
				if(rs[i] == NULL)
				{
					rs[i] = new rational_function_legendre_1d(dimX(), np, nq);
					
					vec _min = min();
					vec _max = max();
					for(int k=0; k<dimX(); ++k)
					{
						if(_min[k] == _max[k])
						{
							_min[k] -= 1.0;
							_max[k] += 1.0;
						}
					}

					rs[i]->setMin(_min) ;
					rs[i]->setMax(_max) ;
				}
				return rs[i];
			}
			else
			{
				std::cout << "<<ERROR>> tried to access out of bound 1D RF" 
				          << std::endl;
				return NULL;
			}
		}

		//! Update the y-1D function for the ith dimension.
		//! \note It will test if the 1D function provided is of the dynamic type
		//! \name rational_function_legendre_1d
		virtual void update(int i, rational_function_1d* r)
		{
			if(dynamic_cast<rational_function_legendre_1d*>(r) != NULL)
			{
				rational_function::update(i, r);
			}
			else
			{
#ifdef DEBUG
				std::cerr << "<<ERROR>> the function provided is not of type \"rational_function_legendre\"" << std::endl;
#endif
			}
		}
} ;

