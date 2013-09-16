#pragma once

#include <functional>
#include <string>
#include <fstream>

#include "common.h"
#include "args.h"
#include "params.h"
#include "data.h"

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
class function : public parametrized
{
	public: // methods

		/* INTERFACE */

		// Overload the function operator
		virtual vec operator()(const vec& x) const = 0 ;
		virtual vec value(const vec& x) const = 0 ;

		//! Load function specific files
		virtual void load(std::istream& in) = 0 ;

		//! \brief Provide a first rough fit of the function. 
		//!
		//! \details
		//! Can be used to set the diffuse component of the function for
		//! example.
		virtual void bootstrap(const data* d, const arguments& args) {}

		//! Provide the dimension of the input space of the function
		virtual int dimX() const { return _nX ; }
		//! Provide the dimension of the output space of the function
		virtual int dimY() const { return _nY ; }

		//! Set the dimension of the input space of the function
		virtual void setDimX(int nX) { _nX = nX ; }
		//! Set the dimension of the output space of the function
		virtual void setDimY(int nY) { _nY = nY ; }

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min) 
		{
#ifdef DEBUG
			assert(min.size() == _nX) ;
#endif
			_min = min ; 
		}
		virtual void setMax(const vec& max) 
		{
#ifdef DEBUG
			assert(max.size() == _nX) ;
#endif
			_max = max ; 
		}
		virtual vec getMin() const { return _min ; }
		virtual vec getMax() const { return _max ; }


		/* EXPORT FUNCTIONS */

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
		double L2_distance(const data* d) const ;

		//! \brief Linf norm to data.
		double Linf_distance(const data* d) const ;


	protected: // data

		// Dimension of the function & domain of definition.
		int _nX, _nY ;
		vec _min, _max ;
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
 * of the derivative and a numerical evaluation of it can be provided.
 */
class nonlinear_function: public function
{
	public: // methods

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
		//
		// The x input of this function is the position in the input space and 
		// has size dimX(), the resulting vector has the size of the parameters
		// times the size of the output domain.
		//
		// The result vector should be orderer as res[i + dimY()*j], output
		// dimension first, then parameters.
		virtual vec parametersJacobian(const vec& x) const = 0;

		//! \brief default non_linear import. Parse the parameters in order.
		virtual void load(std::istream& in)
		{
			// Parse line until the next comment
			while(in.peek() != '#')
			{
				char line[256];
				in.getline(line, 256);
			}

			// Checking for the comment line #FUNC nonlinear_function_phong
			std::string token;
			in >> token;
			if(token.compare("#FUNC") != 0) 
			{ 
				std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
#ifdef DEBUG
				std::cout << "<<DEBUG>> got: \"" << token << "\"" << std::endl;
#endif
			}

			in >> token;
			if(token.compare("nonlinear_function") != 0)
			{
				std::cerr << "<<ERROR>> parsing the stream. A function name is defined." << std::endl;
				std::cerr << "<<ERROR>> did you forget to specify the plugin used to export?" << std::endl;
			}

			int nb_params = nbParameters();
			vec p(nb_params);
			for(int i=0; i<nb_params; ++i)
			{
				in >> token >> p[i];
			}

			setParameters(p);
		}

		//! \brief default non_linear export. It will dump the parameters in order
		//! but won't assign names for the function nor parameters.
		virtual void save_call(std::ostream& out, const arguments& args) const
		{
			if(!args.is_defined("export"))
			{
				// Dump a #FUNC nonlinear
				out << "#FUNC nonlinear_function" << std::endl;

				// Dump the parameters in order
				vec p = parameters();
				for(int i=0; i<p.size(); ++i)
				{
					out << "param_" << i+1 << "\t" << p[i] << std::endl;
				}
				out << std::endl;
			}

			function::save_call(out, args);
		}
};


class compound_function: public nonlinear_function
{
	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const
		{
			return value(x);
		}
		virtual vec value(const vec& x) const
		{
			vec res(_nY);
			res = vec::Zero(_nY);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				if(fs[i]->input_parametrization() != input_parametrization())
				{
					vec temp_x(fs[i]->dimX());
					params::convert(&x[0], input_parametrization(), fs[i]->input_parametrization(), &temp_x[0]);
					res = res + fs[i]->value(temp_x);
				}
				else
				{
					res = res + fs[i]->value(x);
				}
			}
			return res;
		}

		//! Provide a vector like interface
		virtual void push_back(nonlinear_function* f)
		{
			// Update the input param
			if(input_parametrization() == params::UNKNOWN_INPUT)
			{
				setParametrization(f->input_parametrization());
			}
			
			// Update the output param
			if(output_parametrization() == params::UNKNOWN_OUTPUT)
			{
				setParametrization(f->output_parametrization());
			}
			else if(output_parametrization() != f->output_parametrization())
			{
				std::cerr << "Creating a compound function with different output dimensions, this is not allowed" << std::endl;
				throw;
			}
			
			fs.push_back(f);
		}

		//! \brief Access to the i-th function of the compound
		nonlinear_function* operator[](int i) const
		{
#ifdef DEBUG
			assert(i >= 0 && i < fs.size());
#endif
			return fs[i];
		}

		//! \brief Access to the number of elements in the compound object.
		unsigned int size() const
		{
			return fs.size();
		}

		//! Load function specific files
		virtual void load(std::istream& in)
		{
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->load(in);
			}
		}

		//! \brief Provide a first rough fit of the function. 
		virtual void bootstrap(const ::data* d, const arguments& args) 
		{
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->bootstrap(d, args);
			}
		}

		//! Set the dimension of the input space of the function
		virtual void setDimX(int nX) 
		{
			function::setDimX(nX);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->setDimX(nX);
			}
		}
		//! Set the dimension of the output space of the function
		virtual void setDimY(int nY)
		{
			function::setDimY(nY);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->setDimY(nY);
			}
		}

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min) 
		{
			function::setMin(min);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->setMin(min);
			}
		}
		virtual void setMax(const vec& max) 
		{
			function::setMax(max);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->setMax(max);
			}
		}

		//! Number of parameters to this non-linear function
		virtual int nbParameters() const
		{
			int nb_params = 0;
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				nb_params += fs[i]->nbParameters();
			}
			return nb_params;
		}

		//! Get the vector of parameters for the function
		virtual vec parameters() const
		{
			vec params(nbParameters());
			int current_i = 0;
			for(unsigned int f=0; f<fs.size(); ++f)
			{
				int f_size = fs[f]->nbParameters();

				// Handle when there is no parameters to include
				if(f_size > 0)
				{
					vec f_params = fs[f]->parameters();
					for(int i=0; i<f_size; ++i)
					{
						params[i + current_i] = f_params[i];
					}

					current_i += f_size;
				}
			}

			return params;
		}
		
		//! Get the vector of min parameters for the function
		virtual vec getParametersMin() const
		{
			vec params(nbParameters());
			int current_i = 0;
			for(unsigned int f=0; f<fs.size(); ++f)
			{
				int f_size = fs[f]->nbParameters();

				// Handle when there is no parameters to include
				if(f_size > 0)
				{
					vec f_params = fs[f]->getParametersMin();
					for(int i=0; i<f_size; ++i)
					{
						params[i + current_i] = f_params[i];
					}

					current_i += f_size;
				}
			}

			return params;
		}
		
		//! Get the vector of min parameters for the function
		virtual vec getParametersMax() const
		{
			vec params(nbParameters());
			int current_i = 0;
			for(unsigned int f=0; f<fs.size(); ++f)
			{
				int f_size = fs[f]->nbParameters();

				// Handle when there is no parameters to include
				if(f_size > 0)
				{
					vec f_params = fs[f]->getParametersMax();
					for(int i=0; i<f_size; ++i)
					{
						params[i + current_i] = f_params[i];
					}

					current_i += f_size;
				}
			}

			return params;
		}

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p) 
		{
			int current_i = 0;
			for(unsigned int f=0; f<fs.size(); ++f)
			{
				int f_size = fs[f]->nbParameters();

				// Handle when there is no parameters to include
				if(f_size > 0)
				{
					vec f_params(f_size);
					for(int i=0; i<f_params.size(); ++i)
					{
						f_params[i] = p[i + current_i];
					}

					fs[f]->setParameters(f_params);
					current_i += f_size;
				}
			}
		}

		//! \brief Obtain the derivatives of the function with respect to the 
		//! parameters. 
		//
		// The x input of this function is the position in the input space and 
		// has size dimX(), the resulting vector has the size of the parameters
		// times the size of the output domain.
		//
		// The result vector should be orderer as res[i + dimY()*j], output
		// dimension first, then parameters.
		virtual vec parametersJacobian(const vec& x) const
		{
			int nb_params = nbParameters();
			vec jac(nb_params*_nY);
			jac = vec::Zero(nb_params*_nY);

			int start_i = 0;

			// Export the sub-Jacobian for each function
			for(unsigned int f=0; f<fs.size(); ++f)
			{
				nonlinear_function* func = fs[f];
				int nb_f_params = func->nbParameters(); 

				// Only export Jacobian if there are non-linear parameters
				if(nb_f_params > 0)
				{

					vec temp_x(func->dimX());
					params::convert(&x[0], input_parametrization(), func->input_parametrization(), &temp_x[0]);
					vec func_jac = func->parametersJacobian(temp_x);

					for(int i=0; i<nb_f_params; ++i)
					{
						for(int y=0; y<_nY; ++y)
						{
							jac[y + _nY*(i+start_i)] = func_jac[y + _nY*i];
						}
					}
				}

				start_i += nb_f_params;
			}

			return jac;
		}

		//! \brief can set the input parametrization of a non-parametrized
		//! object. Print an error if it is already defined.
		virtual void setParametrization(params::input new_param)
		{
			if(new_param == params::UNKNOWN_INPUT)
				return;

			parametrized::setParametrization(new_param);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				if(fs[i]->input_parametrization() == params::UNKNOWN_INPUT)
				{
					fs[i]->setParametrization(new_param);
				}
			}
		}

		//! \brief can set the output parametrization of a non-parametrized
		//! function. Throw an exception if it tries to erase a previously
		//! defined one.
		virtual void setParametrization(params::output new_param)
		{
			parametrized::setParametrization(new_param);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->setParametrization(new_param);
			}
		}

		//! \brief save function specific data. This has no use for ALTA export
		//! but allows to factorize the code in the C++ or matlab export by
		//! defining function calls that are common to all the plugins.
		virtual void save_body(std::ostream& out, const arguments& args) const
		{
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				fs[i]->save_body(out, args);
				out << std::endl;
			}

			function::save_body(out, args);
		}

		//! \brief save object specific information. For an ALTA export the
		//! coefficients will be exported. For a C++ or matlab export, the call
		//! to the associated function will be done.
		virtual void save_call(std::ostream& out, const arguments& args) const
		{
			bool is_cpp    = args["export"] == "C++";
			bool is_shader = args["export"] == "shader" || args["export"] == "explorer";
			bool is_matlab = args["export"] == "matlab";

			// This part is export specific. For ALTA, the coefficients are just
			// dumped as is with a #FUNC {plugin_name} header.
			//
			// For C++ export, the function call should be done before hand and
			// the line should look like:
			//   res += call_i(x);
			for(unsigned int i=0; i<fs.size(); ++i)
			{
				if(i != 0 && (is_cpp || is_matlab || is_shader))
				{
					out << "\tres += ";
				}

				fs[i]->save_call(out, args);

				if(is_cpp || is_matlab || is_shader)
				{
					out << ";" << std::endl;
				}
			}

			function::save_call(out, args);
		}

	protected:
		std::vector<nonlinear_function*> fs;

};

/*! \brief A Fresnel interface
 *  \ingroup core
 */
class fresnel : public nonlinear_function
{
	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const
		{
			return value(x);
		}
		virtual vec value(const vec& x) const
		{
			vec fres = fresnelValue(x);
			vec func = f->value(x);

			return product(fres, func);
		}

		//! Load function specific files
		virtual void load(std::istream& in)
		{
			if(f != NULL)
			{
				f->load(in);
			}
			else
			{
				std::cout << "<<ERROR>> trying to load a Fresnel object with no base class" << std::endl;
			}
		}

		//! \brief Provide a first rough fit of the function. 
		virtual void bootstrap(const data* d, const arguments& args) 
		{
			fresnelBootstrap(d, args);
			f->bootstrap(d, args);
		}

		//! Set the dimension of the input space of the function
		virtual void setDimX(int nX) 
		{
			function::setDimX(nX);
			f->setDimX(nX);
		}
		//! Set the dimension of the output space of the function
		virtual void setDimY(int nY)
		{
			function::setDimY(nY);
			f->setDimY(nY);
		}

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min) 
		{
			function::setMin(min);
			f->setMin(min);
		}
		virtual void setMax(const vec& max) 
		{
			function::setMax(max);
			f->setMax(max);
		}

		//! Number of parameters to this non-linear function
		virtual int nbParameters() const
		{
			return f->nbParameters() + nbFresnelParameters();
		}

		//! Get the vector of parameters for the function
		virtual vec parameters() const
		{
			int nb_func_params = f->nbParameters();
			int nb_fres_params = nbFresnelParameters();
			int nb_params = nb_func_params + nb_fres_params;

			vec params(nb_params);

			vec func_params = f->parameters();
			for(int i=0; i<nb_func_params; ++i)
			{
				params[i] = func_params[i];
			}

			vec fres_params = getFresnelParameters();
			for(int i=0; i<nb_fres_params; ++i)
			{
				params[i+nb_func_params] = fres_params[i];
			}

			return params;
		}
		
		//! Get the vector of min parameters for the function
		virtual vec getParametersMax() const
		{
			int nb_func_params = f->nbParameters();
			int nb_fres_params = nbFresnelParameters();
			int nb_params = nb_func_params + nb_fres_params;

			vec params(nb_params);

			vec func_params = f->getParametersMax();
			for(int i=0; i<nb_func_params; ++i)
			{
				params[i] = func_params[i];
			}

			vec fres_params = getFresnelParametersMax();
			for(int i=0; i<nb_fres_params; ++i)
			{
				params[i+nb_func_params] = fres_params[i];
			}

			return params;
		}
		
		//! Get the vector of min parameters for the function
		virtual vec getParametersMin() const
		{
			int nb_func_params = f->nbParameters();
			int nb_fres_params = nbFresnelParameters();
			int nb_params = nb_func_params + nb_fres_params;

			vec params(nb_params);

			vec func_params = f->getParametersMin();
			for(int i=0; i<nb_func_params; ++i)
			{
				params[i] = func_params[i];
			}

			vec fres_params = getFresnelParametersMin();
			for(int i=0; i<nb_fres_params; ++i)
			{
				params[i+nb_func_params] = fres_params[i];
			}

			return params;
		}

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p)
		{
			int nb_func_params = f->nbParameters();
			int nb_fres_params = nbFresnelParameters();

			vec func_params(nb_func_params);
			for(int i=0; i<nb_func_params; ++i)
			{
				func_params[i] = p[i];
			}
			f->setParameters(func_params);

			vec fres_params(nb_fres_params);
			for(int i=0; i<nb_fres_params; ++i)
			{
				fres_params[i] = p[i+nb_func_params];
			}
			setFresnelParameters(fres_params);
		}

		//! \brief Obtain the derivatives of the function with respect to the 
		//! parameters.
		virtual vec parametersJacobian(const vec& x) const
		{
			int nb_func_params = f->nbParameters();
			int nb_fres_params = nbFresnelParameters();
			int nb_params = nb_func_params + nb_fres_params;

			vec func_jacobian = f->parametersJacobian(x);
			vec fres_jacobian = getFresnelParametersJacobian(x);

			vec func_value = f->value(x);
			vec fres_value = fresnelValue(x);

			// F = fresnel; f = function
			// d(F * f)(x) /dp = F(x) df(x) /dp + f(x) dF(x) / dp
			vec jac(nb_params*_nY);
			for(int y=0; y<_nY; ++y)
			{
				for(int i=0; i<nb_func_params; ++i)
				{
					jac[y + _nY*i] = func_jacobian[y + _nY*i] * fres_value[y];
				}

				for(int i=0; i<nb_fres_params; ++i)
				{
					jac[y + _nY*(i+nb_func_params)] = fres_jacobian[y + _nY*i] * func_value[y];
				}
			}

			return jac;
		}

		//! \brief set the value for the base function
		void setBase(nonlinear_function* fin)
		{
			f = fin;
		}

		//! \brief provide the input parametrization of the object.
		virtual params::input input_parametrization() const
		{
			return f->input_parametrization();
		}

		//! \brief provide the outout parametrization of the object.
		virtual params::output output_parametrization() const
		{
			return f->output_parametrization();
		}

		//! \brief can set the input parametrization of a non-parametrized
		//! object. Print an error if it is already defined.
		virtual void setParametrization(params::input new_param)
		{
			function::setParametrization(new_param);
			f->setParametrization(new_param);
		}

		//! \brief can set the output parametrization of a non-parametrized
		//! function. Throw an exception if it tries to erase a previously
		//! defined one.
		virtual void setParametrization(params::output new_param)
		{
			function::setParametrization(new_param);
			f->setParametrization(new_param);
		}

	protected: // methods

		//! \brief the interface for the Fresnel code
		virtual vec fresnelValue(const vec& x) const  = 0;

		//! Number of parameters to this non-linear function
		virtual int nbFresnelParameters() const = 0;

		//! Get the vector of parameters for the function
		virtual vec getFresnelParameters() const = 0;

		//! Get the vector of min parameters for the function
		virtual vec getFresnelParametersMin() const
		{
			vec m(nbFresnelParameters());
			for(int i=0; i<nbFresnelParameters(); ++i)
			{
				m[i] = -std::numeric_limits<double>::max();
			}
			return m;
		}

		//! Get the vector of min parameters for the function
		virtual vec getFresnelParametersMax() const
		{
			vec M(nbFresnelParameters());
			for(int i=0; i<nbFresnelParameters(); ++i)
			{
				M[i] = std::numeric_limits<double>::max();
			}
			return M;
		}

		//! Update the vector of parameters for the function
		virtual void setFresnelParameters(const vec& p) = 0;

		//! \brief Obtain the derivatives of the function with respect to the 
		//! parameters.
		virtual vec getFresnelParametersJacobian(const vec& x) const = 0;		

		//! \brief Boostrap the function by defining the diffuse term
		virtual void fresnelBootstrap(const data* d, const arguments& args) = 0;

	protected: //data

		//! the base object
		nonlinear_function* f;
};
