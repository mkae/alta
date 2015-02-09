/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <core/data.h>
#include <core/common.h>
#include <core/args.h>

class data_merl : public data
{
	public: // methods

        data_merl();

		// Load data from a file
		virtual void load(const std::string& filename) ;
		virtual void load(const std::string& filename, const arguments& args) ;
		
		virtual void save(const std::string& filename) const ;

		// Acces to data
		virtual vec get(int i) const;
		virtual vec operator[](int i) const;

      virtual vec value(const vec& in) const;
		
		// Set data
		virtual void set(int i, const vec& x);
		virtual void set(const vec& x);

		// Get data size, e.g. the number of samples to fit
		virtual int size() const ;

		// Get min and max input space values
		virtual vec min() const ;
		virtual vec max() const ;

		virtual int dimX() const ; 
		virtual int dimY() const ; 

	private: // data
		double *brdf ;

} ;
