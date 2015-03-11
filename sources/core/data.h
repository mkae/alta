/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <string>
#include <utility>
#include <iostream>
#include <limits>
#include <fstream>

#include "common.h"
#include "args.h"
#include "params.h"
#include "clustering.h"
#include "ptr.h"

/*! \brief A data object. Allows to load data from files.
 *  \ingroup core
 */
class data : public parametrized
{
	public: // methods

		// Virtual destructor
		virtual ~data() {} 

		// Load data from a file
		virtual void load(const std::string& filename) = 0 ;
		virtual void load(const std::string& filename, const arguments& args) = 0 ;

		// Save the data to a file
		virtual void save(const std::string& filename) const
		{
			std::ofstream file(filename.c_str(), std::ios_base::trunc);
			file << "#DIM " << _nX << " " << _nY << std::endl;
			file << "#PARAM_IN  " << params::get_name(input_parametrization())  << std::endl;
			file << "#PARAM_OUT " << params::get_name(output_parametrization()) << std::endl;
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

    //! \brief Provide an evaluation of the data using interpolation. If
    //! the data object does not provide an interpolation mechanism, it
    //! should throw an exception.
    //!
    //! \details
    //! The input vector must have the parametrization of the data, and
    //! match the total dimension: dimX + dimY.
    virtual vec value(const vec& in) const = 0;

		//! \brief Put the sample inside the data
		virtual void set(const vec& x) = 0;
		virtual void set(int i, const vec& x) = 0;


		// Get data size, e.g. the number of samples to fit
		virtual int size() const = 0 ;

	protected: // data
} ;

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
		data_params(const ptr<data> d, params::input new_param,
		            data_params::clustrering method = data_params::NONE) :
			_clustering_method(method)
		{
			setParametrization(new_param);
			setParametrization(d->output_parametrization());

			_nX = params::dimension(new_param);
			_nY = d->dimY();

			std::cout << "<<INFO>> Reparametrization of the data" << std::endl;
			//TODO
			//clustering<data>(d, _nY, d->parametrization(), new_param, _data);

			std::cout << "<<INFO>> clustering left " << _data.size() << "/" << d->size() << " elements" << std::endl;
			save(std::string("cluster.gnuplot"));
		}

		virtual vec value(const vec&) const
		{
			NOT_IMPLEMENTED();
		}

		// Load data from a file
		virtual void load(const std::string&)
		{
			std::cerr << "<<ERROR>> this data type cannot load data" << std::endl;
			throw;
		}

		virtual void load(const std::string&, const arguments&)
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
		virtual void set(const vec& x)
		{
			this->set(x);
		}

		virtual void set(int i, const vec& x)
		{
			this->set(i, x);
		}

		// Get data size, e.g. the number of samples to fit
		virtual int size() const
		{
			return _data.size();
		}

	protected: // data

		data_params::clustrering _clustering_method;

		std::vector<vec> _data;
};
