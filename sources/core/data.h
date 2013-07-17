#pragma once

#include <string>
#include <utility>
#include <iostream>
#include <limits>
#include <fstream>

#include <QtPlugin>

#include "common.h"
#include "args.h"
#include "params.h"
#include "clustering.h"

/*! \brief A data object. Allows to load data from files.
 *  \ingroup core
 */
class data : public parametrized
{
	public: // methods
/*
		data() : _in_param(params::UNKNOWN_INPUT), _out_param(params::UNKNOWN_OUTPUT) { }
*/

		// Load data from a file
		virtual void load(const std::string& filename) = 0 ;
		virtual void load(const std::string& filename, const arguments& args) = 0 ;

        // Save the data to a file
        virtual void save(const std::string& filename) const
        {
            std::ofstream file(filename.c_str(), std::ios_base::trunc);
            file << "#DIM " << _nX << " " << _nY << std::endl;
            for(int i=0; i<size(); ++i)
            {
                vec x = this->get(i);
                for(int j=0; j<_nX+_nY; ++j)
                {
                    file << x[j] << "\t";
                }
                file << std::endl;
            }
            file.close();
        }

		// Acces to data
		virtual vec get(int i) const = 0 ;
		virtual vec operator[](int i) const = 0 ;

		//! \brief Provide an evaluation in a BRDF way of the data. 
		//!
		//! \details
		//! The input vector in (can have texture coordinate) and the output
		//! vector out are taken to grad a value and return it. The two vectors 
		//! should be compliant with the size and parametrization of the data.
		virtual vec value(vec in, vec out) const = 0;

		//! \brief Put the sample inside the data
		virtual void set(vec x) = 0;


		// Get data size, e.g. the number of samples to fit
		virtual int size() const = 0 ;


		// Get min and max input space values
		virtual vec min() const = 0 ;
		virtual vec max() const = 0 ;


		// Get the size of the input X domain and output Y domain
		virtual int dimX() const { return _nX ; }
		virtual int dimY() const { return _nY ; }
/*
		//! \brief provide the parametrization of the function.
		//! \note some function type can modify the parametrization to adapt
		//! to the data.
		virtual params::input parametrization() const
		{
			return _in_param;
		}

		//! \brief can set the input parametrization of a non-parametrized
		//! function. Throw an exception if it tries to erase a previously
		//! defined one.
		virtual void setParametrization(params::input new_param)
		{
			if(_in_param != params::UNKNOWN_INPUT)
				throw("A parametrization is already defined");

			_in_param = new_param;
		}
*/
	protected: // data

		// Dimensions of the data
		int _nX, _nY ;
/*
		// Input and output parametrization
		params::input  _in_param ;
		params::output _out_param ;
		*/
} ;
	 
Q_DECLARE_INTERFACE(data, "Fitter.Data")

/*! \brief Change the parametrization of data to fit the parametrization of the
 *  function to be fitted.
 *
 *  \ingroup core
 *  \internal
 *  \todo Finish this class
 */
class data_params : public data
{
	public: // structures

		//! \brief when changing from a parametrization to another, you might
		//! lose some dimensions. This list enumerate the different operators
		//! that can be applied on the raw data to be clusterized.
		//! \note by default we use <em>none</em>, but if the input space
		//! dimension is reduced, the program will halt.
		enum clustrering
		{
			MEAN,
			MEDIAN,
			NONE
		};

	public: // methods

		//! \brief contructor requires the definition of a base class that
		//! has a parametrization, and a new parametrization.
		data_params(const data* d, params::input new_param,
				data_params::clustrering method = data_params::NONE) :
			_clustering_method(method)
		{
			_nX = params::dimension(new_param);
			_nY = d->dimY();

			std::cout << "<<INFO>> Reparametrization of the data" << std::endl;
			clustering<data>(d, _nY, d->parametrization(), new_param, _data);

			std::cout << "<<INFO>> clustering left " << _data.size() << "/" << d->size() << " elements" << std::endl;
			save(std::string("cluster.gnuplot"));
		}

		virtual vec value(vec in, vec out) const 
		{
			NOT_IMPLEMENTED();
		}

		// Load data from a file
		virtual void load(const std::string& filename)
		{
			std::cerr << "<<ERROR>> this data type cannot load data" << std::endl;
			throw;
		}

		virtual void load(const std::string& filename, const arguments& args)
		{
			std::cerr << "<<ERROR>> this data type cannot load data" << std::endl;
			throw;
		}

		// Acces to data
		virtual vec get(int i) const
		{
            return _data[i];
		}
		virtual vec operator[](int i) const
		{
			return this->get(i);
		}
		
		//! \todo This should crash at execution.
		virtual void set(vec x)
		{
			this->set(x);
		}

		// Get data size, e.g. the number of samples to fit
		virtual int size() const
		{
            return _data.size();
		}

		// Get min and max input space values
		virtual vec min() const
		{
            return _min;
		}
		virtual vec max() const
		{
            return _max;
		}

	protected: // data

		data_params::clustrering _clustering_method;

        std::vector<vec> _data;

        vec _min, _max;
};
