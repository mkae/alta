#pragma once

#include "function.h"
#include "data.h"
#include "args.h"
#include "common.h"

#include <QtPlugin>

/*! \brief Fitting interface for generic fitting algorithms
 *  \ingroup core
 *
 */
class fitter
{
	public:
		
        //! \brief static function to fit a data set d with the underling
        //! function class. Return the best fit (along with fitting
        //! information ?)
        virtual bool fit_data(const data* d, function* f, const arguments& args) = 0 ;


		virtual void set_parameters(const arguments& args) = 0 ;

} ;

Q_DECLARE_INTERFACE(fitter, "Fitter.Fitter") 
