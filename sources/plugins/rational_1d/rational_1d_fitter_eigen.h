#pragma once

#include "rational_1d_fitter.h"

class rational_1d_fitter_eigen : public rational_1d_fitter
{
	public: // methods

		// Fitting a data object using np elements
		// in the numerator and nq elements in the
		// denominator
		virtual bool fit_data(const data* d, int np, int nq, function*& f) ;

	private:
} ;
