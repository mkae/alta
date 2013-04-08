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
 */
class data_params : public data
{
public:
    // Load data from a file
    virtual void load(const std::string& filename)
    {
        _d->load(filename);
    }

    virtual void load(const std::string& filename, const arguments& args)
    {
        _d->load(filename, args);
    }

    // Acces to data
    virtual vec get(int i) const
    {
        return _d->get(i);
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
    data* _d;
};
