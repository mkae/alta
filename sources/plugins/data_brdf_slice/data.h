#pragma once

#include <core/data.h>
#include <core/common.h>
#include <core/args.h>

class data_brdf_slice : public data
{
	public: // methods

		data_brdf_slice();
		~data_brdf_slice();

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

		// Get min and max input space values
		virtual vec min() const ;
		virtual vec max() const ;

		virtual int dimX() const ; 
		virtual int dimY() const ; 

	private: // data
		double* _data;
		int width, height;
		int slice;

} ;
