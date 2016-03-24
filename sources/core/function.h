/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013 CNRS
   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

//#include <functional>
#include <string>
#include <fstream>
#include <utility>

#include "common.h"
#include "args.h"
#include "params.h"
#include "data.h"

namespace alta {

/*! \brief A representation of an analytical function.
 *  \ingroup core
 *
 *  \details
 *  function are functors with a domain of definition specified by a vector 
 *  interval \f$[\vec{min} .. \vec{max}]\f$ where \f$\vec{min}\f$ and 
 *  \f$\vec{max}\f$ have the size of the input domain.
 *
 *  Any function used by the fitting algorithm should overload publicly this
 *  interface.
 */
class function
{
	public: // methods

    function()
        : _in_param(params::UNKNOWN_INPUT), _out_param(params::UNKNOWN_OUTPUT),
        _nX(0), _nY(0),
        _min(vec::Zero(0)), _max(vec::Zero(0)) { };


		/* NEEDED FUNCTIONS */

		//! \brief Destructor function needed when using shared_ptr
		virtual ~function() {}

		/* INTERFACE */

		// Overload the function operator
		virtual vec operator()(const vec& x) const { return this->value(x); } ;
		virtual vec value(const vec& x) const = 0 ;

		//! \brief Provide a first rough fit of the function. 
		//!
		//! \details
		//! Can be used to set the diffuse component of the function for
		//! example. The default behavior is to load a function file.
		virtual void bootstrap(const ptr<data>, const arguments& args);

		/* IMPORT/EXPORT FUNCTIONS */

		//! Load function specific files
		virtual bool load(std::istream& in) = 0 ;

		//! \brief Save the current function to a specific file type, args 
		//! can be used to differenciate the type of export.
		//! \see rational_function.cpp for an example
		virtual void save(const std::string& filename, const arguments& args) const;

		//! \brief save the header of the output function file. The header 
		//! should store general information about the fit such as the 
		//! command line used the dimension of the fit. L2 and L_inf distance
		//! could be added here.
		virtual void save_header(std::ostream& out, const arguments& args) const ;

		//! \brief save function specific data. This has no use for ALTA 
		//! exportbut allows to factorize the code in the C++ or matlab 
		//! export by defining function calls that are common to all the 
		//! plugins.
		virtual void save_body(std::ostream& out, const arguments& args) const ;

		//! \brief save object specific information. For an ALTA export the
		//! coefficients will be exported. For a C++ or matlab export, the 
		//! call to the associated function will be done.
		virtual void save_call(std::ostream& out, const arguments& args) const ;


		/* METRIC FUNCTIONS */

		//! \brief L2 norm to data.
		//! \note This distance is only valid with respect to the data sampling.
		//! If the measurement points are not uniformly distributed, then the
		//! metric does not represent the real difference integrated over the
		//! hemisphere.
		double L2_distance(const ptr<data>& d) const ;

		//! \brief Linf norm to data.
		//! \note This distance is only valid with respect to the data sampling.
		//! If the measurement points are not uniformly distributed, then the
		//! metric does not represent the real difference integrated over the
		//! hemisphere.
		double Linf_distance(const ptr<data>& d) const ;


    // Definition domain of the function.
    // TODO: Move this to a 'parameterization' class.
    virtual void setMin(const vec& min) { _min = min; }
    virtual void setMax(const vec& max) { _max = max; }
    int dimX() const { return _nX; }
    int dimY() const { return _nX; }
    virtual void setDimX(int x) { _nX = x; _min.resize(x); _max.resize(x); }
    virtual void setDimY(int y) { _nY = y; }
    void setParametrization(params::input p) { _in_param = p; }
    void setParametrization(params::output p) { _out_param = p; }
    params::input parametrization() const { return _in_param; }
    params::input input_parametrization() const { return _in_param; }
    params::output output_parametrization() const { return _out_param; }
    virtual vec min() const { return _min; }
    virtual vec max() const { return _max; }


protected:
    // Input and output parametrization
    // TODO: Move this to a 'parameterization' class.
    params::input  _in_param;
    params::output _out_param;
    int _nX, _nY;
    vec _min, _max;
};

/*! \brief Non-linear function interface
 *  \ingroup core
 *
 * \details
 * Provide a way to obtain the dérivative of the function with respect to its
 * parameters. If the function \f$f(\vec{x})\f$ is defined for a vector of
 * parameters \f$\vec{a}\f$, the resulting vector is \f$df_i = {df \over 
 * da_i}\f$. 
 *
 * \note It is not necessary to have an analytical formulation
 * of the derivative but at least a numerical evaluation of it has to be
 * provided.
 */
class nonlinear_function: public function
{
	public: // methods

		//! \brief Provide a first rough fit of the function.
		//!
		//! \details
		//! The nonllinear_function overload the function bootstrap call to
		//! allows for loading the parameters vector directly from the command
		//! line: --bootstrap [p1, p2, ..., pn]
		virtual void bootstrap(const ptr<data>, const arguments& args);

		//! Number of parameters to this non-linear function
		virtual int nbParameters() const = 0;

		//! Get the vector of parameters for the function
		virtual vec parameters() const = 0;

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p) = 0;

		//! \brief get the maximum value for all the parameters in a vector
		//! format. The maximum vector and the parameter vector have the same
		//! indexing.
		virtual vec getParametersMax() const;

		//! \brief get the minimum value for all the parameters in a vector
		//! format. The minimum vector and the parameter vector have the same
		//! indexing.
		virtual vec getParametersMin() const;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		//!
		//! \details
		//! The \a x input of this function is the position in the 
		//! input space and  has size dimX(), the resulting vector has the 
		//! size of the parameters times the size of the output domain.
		//!
		//! The result vector should be orderer as res[i + dimY()*j], output
		//! dimension first, then parameters.
		virtual vec parametersJacobian(const vec& x) const = 0;

		//! \brief default non_linear import. Parse the parameters in order.
		virtual bool load(std::istream& in);

		//! \brief default non_linear export. It will dump the parameters in 
		//! order but won't assign names for the function nor parameters.
		virtual void save_call(std::ostream& out, const arguments& args) const;
};

/*! \brief A compound (sum) of multiple nonlinear functions. This allows to
 *  perform the fit of multiple lobes at once.
 *  \ingroup core
 */
class compound_function: public nonlinear_function
{
	public: // methods

		//Destructor
		virtual ~compound_function();

		// Overload the function operator
		virtual vec operator()(const vec& x) const;
		virtual vec value(const vec& x) const;

		//! Provide a vector like interface
		//! This function allows to put a new nonlinear function \a f in the 
		//! compound object. This function will be processed for nonlinear
		//! optimisation only if \a fixed equals true.
		virtual void push_back(const ptr<nonlinear_function>& f, const arguments& f_args);

		//! \brief Access to the i-th function of the compound
		nonlinear_function* operator[](int i) const;

		//! \brief Access to the number of elements in the compound object.
		unsigned int size() const;

		//! Load function specific files
		virtual bool load(std::istream& in);

		//! \brief Provide a first rough fit of the function. 
		//! For compound object, you can define the first guess using the
		//! either the global function or using the individual command per
		//! function. <br />
		//!
		//! <u>Examples:</u><br />
		//! \verbatim
		//! --func [libfunc1.so, libfunc2.so] --bootstrap first_guess.brdf
		//! \endverbatim
		//! Will load the file <em>first_guess.brdf</em> as the initial value
		//! <br />
		//! \verbatim
		//! --func [libfunc1.so --bootstrap [val1, val2], libfunc2.so --bootstrap first_guess1.brdf]
		//! \endverbatim
		//! Will load the vector of parameters <em>[val1, val2]</em> for the
		//! first function and the file <em>first_guess1.brdf</em> for the
		//! second one. <br />
		//!
		//! <u>Local/Global policy:</u></br />
		//! Local bootstrap can not overload the global bootstrap.
		virtual void bootstrap(const ptr<data> d, const arguments& args);

		//! Set the dimension of the input space of the function
		virtual void setDimX(int nX);

		//! Set the dimension of the output space of the function
		virtual void setDimY(int nY);

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min);
		virtual void setMax(const vec& max);

		//! Number of parameters to this non-linear function
		virtual int nbParameters() const;

		//! Get the vector of parameters for the function
		virtual vec parameters() const;
		
		//! Get the vector of min parameters for the function
		virtual vec getParametersMin() const;
		
		//! Get the vector of min parameters for the function
		virtual vec getParametersMax() const;

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p);

		//! \brief Obtain the derivatives of the function with respect to the 
		//! parameters. 
		//
		// The x input of this function is the position in the input space and 
		// has size dimX(), the resulting vector has the size of the parameters
		// times the size of the output domain.
		//
		// The result vector should be orderer as res[i + dimY()*j], output
		// dimension first, then parameters.
		virtual vec parametersJacobian(const vec& x) const;

		//! \brief can set the input parametrization of a non-parametrized
		//! object. Print an error if it is already defined.
		virtual void setParametrization(params::input new_param);

		//! \brief can set the output parametrization of a non-parametrized
		//! function. Throw an exception if it tries to erase a previously
		//! defined one.
		virtual void setParametrization(params::output new_param);

		//! \brief save function specific data. This has no use for ALTA export
		//! but allows to factorize the code in the C++ or matlab export by
		//! defining function calls that are common to all the plugins.
		virtual void save_body(std::ostream& out, const arguments& args) const;

		//! \brief save object specific information. For an ALTA export the
		//! coefficients will be exported. For a C++ or matlab export, the call
		//! to the associated function will be done.
		virtual void save_call(std::ostream& out, const arguments& args) const;

	protected:
		std::vector<ptr<nonlinear_function>> fs;
		std::vector<arguments> fs_args;
		std::vector<bool> is_fixed;

};

/*! \brief A product of nonlinear functions. This class aims to simplify the
 *  fitting and progressive fitting of Fresnel effects that can be separated
 *  from the main lobe.
 *  \ingroup core
 */
class product_function : public nonlinear_function
{
	public: // methods

		//! \brief Constructor of the product function, affect the two function
		//! to already created nonlinear_function objects.
		product_function(const ptr<nonlinear_function>& g1, const ptr<nonlinear_function>& g2, 
                         bool is_g1_fixed = false, bool is_g2_fixed = false);

		~product_function();

		/* ACCESS TO INDIVIDUAL ELEMENTS */

		//! \brief Access to the first member of the product
		nonlinear_function* first() const;

		//! \biref Access to the second member of the product
		nonlinear_function* second() const;


		/* EVALUATION FUNCTIONS */

		//! \brief Overload the function operator, directly call the value function.
		virtual vec operator()(const vec& x) const;

		//! \brief Return the product between the value of function f1 and function
		//! f2. The input parameter x should be in the parametrization of f1. This
		//! function will do the conversion before getting f2's value.
		virtual vec value(const vec& x) const;


		/* IMPORT/EXPORT FUNCTIONS */
		
		//! \brief Load the two functions in order f1, then f2. If one of the 
		//! function cannot be loaded, this function will continue. It will only
		//! return false if both functions cannot be loaded.
		virtual bool load(std::istream& in);

		//! \brief save function specific data. This has no use for ALTA export
		//! but allows to factorize the code in the C++ or matlab export by
		//! defining function calls that are common to all the plugins.
		virtual void save_body(std::ostream& out, const arguments& args) const;

		//! \brief save object specific information. For an ALTA export the
		//! coefficients will be exported. For a C++ or matlab export, the call
		//! to the associated function will be done.
		virtual void save_call(std::ostream& out, const arguments& args) const;
		

		//! \brief Provide a first rough fit of the function. 
		virtual void bootstrap(const ptr<data> d, const arguments& args);

		// Set the dimension of the input/output space of the function
		virtual void setDimX(int nX);
		virtual void setDimY(int nY);

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min);
		virtual void setMax(const vec& max);
		
		//! Provide the output parametrization of the object.
		virtual params::output output_parametrization() const;
		
		//! Set the input/output parametrization of a non-parametrized
		//! object. Print an error if it is already defined.
		virtual void setParametrization(params::input  new_param);
		virtual void setParametrization(params::output new_param);

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const;
		
		//! \biref Get the vector of parameters for the function
		virtual vec parameters() const;

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p);

		// Get the vector of min/max parameters for the function
		virtual vec getParametersMax() const;
		virtual vec getParametersMin() const;

		//! \brief Obtain the derivatives of the function with respect to the 
		//! parameters.
		virtual vec parametersJacobian(const vec& x) const;

	private: // data

		// Function composing the product
		ptr<nonlinear_function> f1, f2;
		
		std::pair<bool,bool> _is_fixed; /*!< represents whether or not the parameters of each function is fixed regardint the optimizer */

};

class cosine_function : public nonlinear_function
{
	public:
		// Set the input parametrization to CARTESIAN to reduce the number
		// of transformations in a compound object.
    cosine_function()
    {
        _in_param = params::CARTESIAN;
        _nX = 6;
    }

		// Overload the function operator
		virtual vec operator()(const vec& x) const 
		{
			return value(x);
		}
		virtual vec value(const vec& x) const
		{
			vec res(dimY());
			for(int i=0; i<dimY(); ++i) { res[i] = ((x[2] > 0.0) ? x[2] : 0.0) * ((x[5] > 0.0) ? x[5] : 0.0); }
			return res;
		}

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const
		{
			return 0;
		}

		//! Get the vector of parameters for the function
		vec parameters() const 
		{
			vec res(1);
			return res;
		}

		//! Update the vector of parameters for the function
		void setParameters(const vec&) 
		{
		}

		//! Obtain the derivatives of the function with respect to the 
		//! parameters. 
		vec parametersJacobian(const vec&) const 
		{
			vec jac(1);
			return jac;
		}
};
}
