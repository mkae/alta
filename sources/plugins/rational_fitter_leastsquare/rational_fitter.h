#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/rational_function.h>
#include <core/vertical_segment.h>

class rational_fitter_leastsquare : public fitter
{
	public: // methods

		rational_fitter_leastsquare() ;
		virtual ~rational_fitter_leastsquare() ;

		// Fitting a data object
		//
		virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& fit) ;
		virtual bool fit_data(const ptr<vertical_segment>& dat, int np, int nq, int ny, rational_function_1d* fit) ;

	protected: // data

		int _np, _nq ;
		int _max_iter;
} ;

