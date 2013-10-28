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
		rational_function_legendre_1d(int np, int nq) ;
		rational_function_legendre_1d(const vec& a, const vec& b) ;
		virtual ~rational_function_legendre_1d() {}

		// Get the p_i and q_j function
		virtual double p(const vec& x, int i) const ;
		virtual double q(const vec& x, int j) const ;

	protected:  // methods

		// Legendre polynomial evaluation
		double legendre(double x, int i) const;
} ;

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
					rs[i] = new rational_function_legendre_1d(np, nq);
					rs[i]->setDimX(dimX());
					rs[i]->setDimY(dimY());

					// Test if the input domain is not empty. If one dimension of
					// the input domain is a point, I manually inflate this dimension
					// to avoid numerical issues.
					vec _min = getMin();
					vec _max = getMax();
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

