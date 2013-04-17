#pragma once

#include <string>
#include <utility>

#include <QtPlugin>

#include "common.h"
#include "args.h"
#include "params.h"

/*! \brief A data object. Allows to load data from files.
 *  \ingroup core
 */
class data
{
		public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) = 0 ;
		virtual void load(const std::string& filename, const arguments& args) = 0 ;

		// Acces to data
		virtual vec get(int i) const = 0 ;
		virtual vec operator[](int i) const = 0 ;
		virtual vec value(vec in, vec out) const 
		{
			return vec(_nY) ;
		}

		// Get data size, e.g. the number of samples to fit
		virtual int size() const = 0 ;

		// Get min and max input space values
		virtual vec min() const = 0 ;
		virtual vec max() const = 0 ;

        virtual int dimX() const { return _nX ; }
        virtual int dimY() const { return _nY ; }

        virtual params::input parametrization() const
        {
            return params::UNKNOWN_INPUT;
        }
	
	protected:
		// Dimension of the function
		int _nX, _nY ;
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
        mean,
        median,
        none
    };

public: // methods

    //! \brief contructor requires the definition of a base class that
    //! has a parametrization, and a new parametrization.
    data_params(const data* d, params::input param,
                data_params::clustrering method = data_params::none) :
        _d(d), _param_in(param), _clustering_method(method)
    {
        _nX = params::dimension(param);
        _nY = d->dimY();

        if(_nX < _d->dimX() && method == data_params::none)
        {
            throw("No cluster method provided");
        }
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
        vec res(_nX + _nY);
        vec in = _d->get(i);

        params::convert(&in[0], _d->parametrization(), _param_in, &res[0]);
        memcpy(&res[_nX], &in[_d->dimX()], _nY*sizeof(double));

        return res;
    }
    virtual vec operator[](int i) const
    {
        return this->get(i);
    }

    // Get data size, e.g. the number of samples to fit
    virtual int size() const
    {
        return _d->size();
    }

    // Get min and max input space values
    virtual vec min() const
    {
        return _d->min();
    }
    virtual vec max() const
    {
        return _d->max();
    }

protected:
    // data object to interface
    const data* _d;
    params::input _param_in;
    data_params::clustrering _clustering_method;
};
