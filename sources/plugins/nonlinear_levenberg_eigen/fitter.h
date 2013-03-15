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
#include <core/vertical_segment.h>

class nonlinear_fitter_eigen: public QObject, public fitter
{
	Q_OBJECT
	Q_INTERFACES(fitter)

	public: // methods
	
		nonlinear_fitter_eigen() ;
		virtual ~nonlinear_fitter_eigen() ;
			
		// Fitting a data object
		//
		virtual bool fit_data(const data* d, function* fit) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

		// Unsable stuff
		virtual data*     provide_data() const ;
		virtual function* provide_function() const ;

	protected: // function


	protected: // data
} ;

