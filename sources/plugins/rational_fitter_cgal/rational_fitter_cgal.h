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

class rational_fitter_cgal : public QObject, public fitter
{
	Q_OBJECT
	Q_INTERFACES(fitter)

	public: // methods
	
		rational_fitter_cgal() ;
		virtual ~rational_fitter_cgal() ;

		// Fitting a data object
		//
		virtual bool fit_data(const data* d, function* fit) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

		// Obtain associated data and functions
		//
		virtual data*     provide_data() const ;
		virtual function* provide_function() const ;

	protected: // data

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const rational_data* d, int np, int nq, rational_function* fit) ;
		virtual bool fit_data(const rational_data* dat, int np, int nq, int ny, rational_function* fit) ;

		// min and Max usable np and nq values for the fitting
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;
} ;
