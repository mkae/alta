#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <QObject>
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>

#include <rational_function.h>
#include <rational_data.h>


class rational_fitter_quadprog : public QObject, public fitter
{
	Q_OBJECT
	Q_INTERFACES(fitter)

	public: // methods
	
		rational_fitter_quadprog() ;
		virtual ~rational_fitter_quadprog() ;

		// Fitting a data object
		virtual bool fit_data(const data* d, function* fit) ;

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const data* d, int np, int nq, function* fit) ;

		virtual void set_parameters(const arguments& args) ;

		virtual data* provide_data() const
		{
			return new rational_data() ;
		}	;
		virtual function* provide_function() const 
		{
			return new rational_function() ;
		}	;
	protected: // data

		// min and Max usable np and nq values for the fitting
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;
} ;

