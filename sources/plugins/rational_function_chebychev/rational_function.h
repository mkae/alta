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

class rational_function_chebychev_1d : public rational_function_1d
{
public: // methods

    rational_function_chebychev_1d() ;
    rational_function_chebychev_1d(int np, int nq) ;
    rational_function_chebychev_1d(const vec& a, const vec& b) ;
    virtual ~rational_function_chebychev_1d() {}

    // Get the p_i and q_j function
    virtual double p(const vec& x, int i) const ;
    virtual double q(const vec& x, int j) const ;

protected:  // methods

} ;

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

