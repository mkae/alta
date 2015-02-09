/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <QObject>
#include <core/data.h>
#include <core/common.h>
#include <core/args.h>

class data_astm : public QObject, public data
{
//	Q_OBJECT
   Q_INTERFACES(data)

	public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) ;
		virtual void load(const std::string& filename, const arguments& args) ;

		// Acces to data
		virtual vec get(int i) const ;
		virtual vec& operator[](int i) ;
		virtual vec value(vec in, vec out) const ;

		// Get data size, e.g. the number of samples to fit
		virtual int size() const ;

		// Get min and max input space values
		virtual vec min() const ;
		virtual vec max() const ;

		virtual int dimX() const
		{
			return 4;
		}
		virtual int dimY() const
		{
			return 3;
		}	

		virtual params::input parametrization() const
		{
			return params::SPHERICAL_TL_PL_TV_PV;
		}


	private: // data
		std::vector<vec> _data;

} ;
