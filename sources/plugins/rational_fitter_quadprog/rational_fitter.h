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

/*! \todo This plugin is not working with 1D example. This is weird. Use the following command
 *  to generate the bug:
     \verbatim 
     ./build/plugin_loader --input ../data/1d/transmi_lame_50_50_avg45.txt --output output.rational --fitter librational_fitter_quadprog.so --min 400 --max 500 --min-np 1 --min-nq 1 --np 20 --nq 20 --dt 0.01
     \endverbatim
     It should freeze after np = 6 and nq = 6
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
		virtual bool fit_data(const vertical_segment* d, int np, int nq, rational_function* fit) ;
		virtual bool fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function* fit) ;

		// min and Max usable np and nq values for the fitting
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;
} ;

