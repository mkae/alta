#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <QObject>
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/fitter.h>
#include <core/args.h>

/*! \brief A rational fitter using the library QuadProg++
 *  You can find the library here: http://quadprog.sourceforge.net/
 */
class rational_fitter_quadprog : public QObject, public fitter
{
	Q_OBJECT
	Q_INTERFACES(fitter)

	public: // methods
	
		rational_fitter_quadprog() ;
		virtual ~rational_fitter_quadprog() ;

		// Fitting a data object
		//
        virtual bool fit_data(const data* d, function* fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

		// Obtain associated data and functions
		//
		virtual data*     provide_data() const ;
		virtual function* provide_function() const ;

    private: // methods

	protected: // data

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const vertical_segment* d, int np, int nq, rational_function* fit) ;
		virtual bool fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function_1d* fit) ;

		// min and Max usable np and nq values for the fitting
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;

        // Add constraints to the boundary of the domain. You can shrink it of
        // the parameter --boundary-constraint *double*
        double _boundary;
} ;

