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

/*! \brief A fitter for rational functions that can run in parallel
 *  using QuadProg++ quadratic solver.
 *  \ingroup plugins
 *
 *  \details
 *  You can find QuadProg++ here: http://quadprog.sourceforge.net/
 */
class rational_fitter_parallel : public QObject, public fitter
{
	Q_OBJECT
	Q_INTERFACES(fitter)

	public: // methods
	
		rational_fitter_parallel() ;
		virtual ~rational_fitter_parallel() ;

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
		virtual bool fit_data(const vertical_segment* d, int np, int nq, rational_function* fit, vec& p, vec& q, double& delta) ;
		virtual bool fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function* fit, vec& p, vec& q, double& delta) ;

        //! \brief Create a constraint vector given its index i in the data
        // object and the rational function object to fit.
        virtual void get_constraint(int i, int np, int nq, int ny, const vertical_segment* data, const rational_function* func, vec& cu, vec& cl);

		// min and Max usable np and nq values for the fitting
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;
} ;

