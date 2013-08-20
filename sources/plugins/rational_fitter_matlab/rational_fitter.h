#pragma once

// Include STL
#include <vector>
#include <string>
#include <engine.h>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/rational_function.h>
#include <core/vertical_segment.h>

class rational_fitter_matlab : public fitter
{
	public: // methods

		rational_fitter_matlab() ;
		virtual ~rational_fitter_matlab() ;

		// Fitting a data object
		//
		virtual bool fit_data(const data* d, function* fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // function

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const vertical_segment* d, int np, int nq, rational_function* fit) ;
		virtual bool fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function_1d* fit) ;

	protected: // data

		// min and Max usable np and nq values for the fitting
		// 
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;
		Engine *ep;
} ;

