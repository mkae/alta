#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/vertical_segment.h>

// CERES include
#include <ceres/ceres.h>

/*! \brief A non-linear fitter using the CERES solver
 *  \ingroup plugins
 *
 *  \details
 *  <h3>Third party requirements</h3>
 *
 *  You will need three external libraries to compile this plugin:
 *  <ul>
 *		<li><a href="https://ceres-solver.googlesource.com/ceres-solver">CERES</a>
 *		library, version 1.5.0</li>
 *		<li><a href="http://code.google.com/p/google-glog">Google glog</a> library
 *		version 0.3.1</li>
 *		<li><a href="http://eigen.tuxfamily.org/">Eigen library</a> version 3</li>
 *  </ul>
 *
 *  The last two dependencies are required to compile CERES and this plugin
 *  must be linked against Google glog to run.
 *
 *  You need to provide your own ceres.prf file for qmake to generate the correct
 *  Makefile or Visual Studio solution. In particular this configuration file
 *  should provide:
 *
 *  <pre>
 *  INCLUDEPATH += [path-to-ceres-include]
 *  LIBS += -L[path-to-ceres-lib] -lceres -L[path-to-glog-lib] -lglog -lgomp
 *  </pre>
 *
 *
 *  <h3>Plugin parameters</h3>
 *
 *  We provide the following command line arguments to manipulate this plugin:
 *  <ul>
 *		<li><b>--ceres-max-num-iterations</b> <em>[int]</em> to control the number
 *		of iterations the non linear solver will take before returning a solution</li>
 *		<li><b>--ceres-factorizer</b> <em>[string]</em> to control the type of dense
 *		factorization method used to solve the <a href="http://homes.cs.washington.edu/~sagarwal/ceres-solver/stable/solving.html?highlight=dense_qr#equation-normal">normal equations</a></li>
 *  </ul>
 */
class nonlinear_fitter_ceres: public fitter
{
	public: // methods

		nonlinear_fitter_ceres() ;
		virtual ~nonlinear_fitter_ceres() ;

		// Fitting a data object
		//
		virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

		// Provide user parameters to the fitter
		//
		virtual void set_parameters(const arguments& args) ;

	protected: // data

        // Fitter options
        ceres::Solver::Options options;
} ;
