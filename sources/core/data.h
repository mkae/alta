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

        virtual params::type parametrization() const
        {
            return params::UNKNOWN;
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
public:
    //! \brief contructor requires the definition of a base class that
    //! has a parametrization, and a new parametrization.
    data_params(const data* d, params::type param) : _d(d), _param(param)
    {
        _nX = d->dimX(); //! \todo the parametrization can change the size
                         //! of the input domain.
        _nY = d->dimY();
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

        params::convert(&_d->get(i)[0], _d->parametrization(), _param, &res[0]);

        return res;
    }
    virtual vec operator[](int i) const
    {
        return (*_d)[i];
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
    params::type _param;
};
