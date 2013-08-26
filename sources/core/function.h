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

		// Overload the function operator
		virtual vec operator()(const vec& x) const = 0 ;
		virtual vec value(const vec& x) const = 0 ;

		//! Load function specific files
		virtual void load(const std::string& filename) = 0 ;

		//! \brief Provide a first rough fit of the function. 
		//!
		//! \details
		//! Can be used to set the diffuse component of the function for
		//! example.
		virtual void bootstrap(const data* d, const arguments& args) {}

		//! \brief Save the current function to a specific file type, args can 
		//! be used to differenciate the type of export.
		//!
		//! \see rational_function.cpp for an example
		virtual void save(const std::string& filename, const arguments& args) const
		{
			std::cout << "<<DEBUG>> Exporting the function" << std::endl;
			if(args.is_defined("export"))
			{
				if(args["export"].compare("c++") == 0)
				{
					std::cout << "<<INFO>> will export in C++ format" << std::endl;
					save_cpp(filename, args);
				}
				else if(args["export"].compare("matlab") == 0)
				{
					std::cout << "<<INFO>> will export in matlab format" << std::endl;
					save_matlab(filename, args);
				}
				else if(args["export"].compare("explorer") == 0)
				{
					std::cout << "<<INFO>> will export in BRDF explorer format" << std::endl;
					save_brdfexplorer(filename, args);
				}
				else
				{
					std::cerr << "<<ERROR>> the export format is unknown" << std::endl ;
				}
			}
			else
			{
				save_alta(filename, args) ;
			}
		}

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

	protected: // function

		//! \brief Standard saving function.
		virtual void save_alta(const std::string& filename, const arguments& args) const 
		{
			NOT_IMPLEMENTED();
		}

		//! \brief Output the function as a gnuplot file. It requires
		//! the data object to output the function at the input location only.
		virtual void save_gnuplot(const std::string& filename, const data* d, 
				const arguments& args) const
		{
#ifndef OLD
			std::ofstream file(filename.c_str(), std::ios_base::trunc);
			for(int i=0; i<d->size(); ++i)
			{
				vec v = d->get(i) ;
				//		vec y1 ; y1.assign(d->dimY(), 0.0) ;
				//		for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

				vec y2 = value(v) ;
				for(int u=0; u<d->dimX(); ++u)
					file << v[u] << "\t" ;

				for(int u=0; u<d->dimY(); ++u)
					file << y2[u] << "\t" ;

				file << std::endl ;
			}
			file.close();
#else
			NOT_IMPLEMENTED();
#endif
		}

		//! \brief Output the function using a C++ function formating.
		virtual void save_cpp(const std::string& filename, const arguments& args) const 

		{
			NOT_IMPLEMENTED();
		}

		//! \brief Output the function using a C++ function formating.
		virtual void save_matlab(const std::string& filename, const arguments& args) const 
		{
			NOT_IMPLEMENTED();
		}

		//! \brief Output the function using a BRDF Explorer formating.
		virtual void save_brdfexplorer(const std::string& filename, const arguments& args) const
		{
			NOT_IMPLEMENTED();
		}

		//! \brief parse the header of the file and return the corresponding
		//! arguments and associate stream
		void load_header(const std::string& filename, arguments& args, std::ifstream& file)
		{
			file.open(filename.c_str()) ;
			if(!file.is_open())
			{
				std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
				throw ;
			}

			while(file.peek() == '#')
			{
				std::string line ;
				std::getline(file, line) ;
				std::stringstream linestream(line) ;

				linestream.ignore(1) ;

				std::string comment ;
				linestream >> comment ;

				if(comment == std::string("DIM"))
				{
					linestream >> _nX >> _nY ;
				}
				else if(comment == std::string("CMD"))
				{
					args = arguments::create_arguments(line.substr(5, std::string::npos));
				}
			}
		}
		
		void save_header(const std::string& filename, arguments& args, std::ofstream& file)
		{
			file.open(filename.c_str()) ;
			if(!file.is_open())
			{
				std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
				throw ;
			}

			file << "#CMD " << args.get_cmd() << std::endl;
			file << "#DIM " << _nX << " " << _nY << std::endl;
			file << "#PARAM_IN  " << params::get_name(input_parametrization()) << std::endl;
			//file << "#PARAM_OUT " << params::get_name(output_parametrization()) << std::endl;
			file << std::endl;
		}

	public: // methods

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
};


class compound_function: public nonlinear_function, public std::vector<nonlinear_function*>
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
			for(int i=0; i<this->size(); ++i)
			{
				res = res + this->at(i)->value(x);
			}
			return res;
		}

		//! Load function specific files
		virtual void load(const std::string& filename)
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->load(filename);
			}
		}
		
		//! \brief Provide a first rough fit of the function. 
		virtual void bootstrap(const ::data* d, const arguments& args) 
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->bootstrap(d, args);
			}
		}
		
		//! Save the Fresnel part along with the function
		virtual void save(const std::string& filename, const arguments& args) const
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->save(filename, args);
			}
		}

		//! Set the dimension of the input space of the function
		virtual void setDimX(int nX) 
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->setDimX(nX);
			}
		}
		//! Set the dimension of the output space of the function
		virtual void setDimY(int nY)
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->setDimY(nY);
			}
		}

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min) 
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->setMin(min);
			}
		}
		virtual void setMax(const vec& max) 
		{
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->setMax(max);
			}
		}

		//! Number of parameters to this non-linear function
		virtual int nbParameters() const
		{
			int nb_params = 0;
			for(int i=0; i<this->size(); ++i)
			{
				nb_params += this->at(i)->nbParameters();
			}
		}

		//! Get the vector of parameters for the function
		virtual vec parameters() const
		{
			vec params(nbParameters());
			int current_i = 0;
			for(int f=0; f<this->size(); ++f)
			{
				vec f_params = this->at(f)->parameters();
				for(int i=0; i<f_params.size(); ++i)
				{
					params[i + current_i] = f_params[i];
				}

				current_i += f_params.size();
			}

			return params;
		}

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p) 
		{
			int current_i = 0;
			for(int f=0; f<this->size(); ++f)
			{
				int f_size = this->at(f)->nbParameters();
				vec f_params(f_size);
				for(int i=0; i<f_params.size(); ++i)
				{
					f_params[i] = p[i + current_i];
				}

				this->at(f)->setParameters(f_params);
				current_i += f_size;
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
			int nb_params = this->nbParameters();
			vec jac(nb_params*_nY);

			jac.assign(nb_params*_nY, 0.0);

			int start_i = 0;

			// Export the sub-Jacobian for each function
			for(int f=0; f<this->size(); ++f)
			{
				nonlinear_function* func = this->at(f);
				int nb_f_params = func->nbParameters(); 

				vec func_jac = func->parametersJacobian(x);

				for(int i=0; i<nb_f_params; ++i)
				{
					for(int y=0; y<_nY; ++y)
					{
						jac[y + _nY*(i+start_i)] = func_jac[y + _nY*i];
					}
				}

				start_i += nb_f_params;
			}
		}
	
		//! \brief can set the input parametrization of a non-parametrized
		//! object. Print an error if it is already defined.
		virtual void setParametrization(params::input new_param)
		{
			parametrized::setParametrization(new_param);
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->setParametrization(new_param);
			}
		}
		
		//! \brief can set the output parametrization of a non-parametrized
		//! function. Throw an exception if it tries to erase a previously
		//! defined one.
		virtual void setParametrization(params::output new_param)
		{
			parametrized::setParametrization(new_param);
			for(int i=0; i<this->size(); ++i)
			{
				this->at(i)->setParametrization(new_param);
			}
		}

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
			return fres * f->value(x);
		}

		//! Load function specific files
		virtual void load(const std::string& filename)
		{
			if(f != NULL)
			{
				f->load(filename);
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
		
		//! Save the Fresnel part along with the function
		virtual void save(const std::string& filename, const arguments& args) const
		{
			f->save(filename, args);
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
			for(int i=nb_func_params; i<nb_params; ++i)
			{
				params[i] = fres_params[i-nb_func_params];
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
					jac[y + _nY*(i+nb_func_params)] = fres_jacobian[y + _nY*i] * fres_value[y];
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
