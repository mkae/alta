#pragma once

#include <core/data.h>
#include <core/common.h>
#include <core/args.h>

#include <flann/flann.hpp>

/*! \brief Load a data file, but provide access to an interpolated version of
 *  the data points.
 *
 *  <h2>Requirements:</h2>
 *  The FLANN library.
 *  On linux plateforms it can be obtained using the package manager: 
 *  \verbatim 
 *		sudo apt-get install libflann-dev
 *	\endverbatim
 */
class data_interpolant : public data
{
	public: // methods

		data_interpolant();
		~data_interpolant();

		// Load data from a file
		virtual void load(const std::string& filename) ;
		virtual void load(const std::string& filename, const arguments& args) ;

		virtual void save(const std::string& filename) const ;

		// Acces to data
		virtual vec get(int i) const ;
		virtual vec operator[](int i) const ;

		virtual vec value(vec in, vec out) const;
		virtual vec value(vec x) const;

		// Set data
		virtual void set(vec x);

		// Get data size, e.g. the number of samples to fit
		virtual int size() const ;

	private: // data
		
		// The data object used to load sparse points sets
		ptr<data> _data;

		// Interpolation 
		flann::Index< flann::L2<double> >* _kdtree;
		flann::Matrix<double>* _points;
		int _knn;
} ;
