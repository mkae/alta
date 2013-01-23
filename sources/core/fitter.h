#pragma once

#include "function.h"
#include "data.h"
#include "args.h"

#include <QtPlugin>

/*
 * Fitting interface for generic fitting algorithms
 *
 */
class fitter
{
	public:
		
		// Static function to fit a data set d with the
		// underling function class. Return the best
		// fit (along with fitting information ?)
		//
		virtual bool fit_data(const data* d, function*& f) = 0 ;


		virtual void set_parameters(const arguments& args) = 0 ;
} ;

Q_DECLARE_INTERFACE(fitter, "Fitter.Fitter") 
