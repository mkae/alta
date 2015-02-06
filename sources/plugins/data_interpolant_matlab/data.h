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

#include <engine.h>


/*! \brief Load a data file, but provide access to an interpolated version of
 *  the data points.
 *
 *  <h2>Requirements:</h2>
 *  Matlab engine
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
		Engine *ep;		
} ;
