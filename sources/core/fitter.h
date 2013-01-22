#pragma once

#include "function.h"
#include "data.h"

/*
 * Fitting interface for generic fitting algorithms
 *
 */
template<class X, class Y> class fitter
{
	public:
		
		// Static function to fit a data set d with the
		// underling function class. Return the best
		// fit (along with fitting information ?)
		//
		virtual bool fit_data(const data<X, Y>& d, function<X, Y>& f) = 0 ;

} ;
