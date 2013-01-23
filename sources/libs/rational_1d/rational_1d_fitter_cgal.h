#pragma once

// Include STL
#include <functional>
#include <vector>
#include <string>
#include <tuple>

// Personal include
#include "rational_1d_fitter.h"


class rational_1d_fitter_cgal : public rational_1d_fitter
{
	public: // methods

		// Fitting a data object
		virtual bool fit_data(const rational_1d_data& data, rational_1d& fit) ;

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator.
		virtual bool fit_data(const rational_1d_data& data, int np, int nq, rational_1d& fit) ;
} ;

