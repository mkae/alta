/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>
#include <core/common.h>
#include <core/args.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

class ASTM : public data
{

private: // data
	std::vector<vec> _data;

public: //methods
	ASTM() : data() {

	}

	// Load data from a file
	virtual void load(const std::string& filename) 
	{
		std::ifstream file(filename.c_str());
		std::string line;
		
		// Parse the header
		for(int i=0; i<22; ++i)
		{
			std::getline(file, line);
		}

		while(file.good())
		{
			std::getline(file, line);
			
			size_t s = 0; // start
			size_t t;     // end
			bool sample_wrong = false;
			vec x(7); 
			for(int i=0; i<7; ++i)
			{
				t = line.find_first_of(',', s);
				if(t == s || t == std::string::npos)
				{
					sample_wrong = true;
				}

				x[i] = atof(line.substr(s, t-s).c_str());

				s = t+1;
			}

			if(!sample_wrong)
			{
				_data.push_back(x);
			}
		}

		file.close();
	}
	virtual void load(const std::string& filename, const arguments& args)
	{
		this->load(filename);
	}

	// Acces to data
	virtual vec get(int i) 
	{
		return _data[i];
	}
	virtual vec operator[](int i) const 
	{
		return _data[i];
	}
	virtual vec value(const vec& in) const
	{
		std::cout << "<<ERROR>> not implemented" << std::endl;
		vec res(4);
		return res;
	}

	// Get data size, e.g. the number of samples to fit
	int size() const 
	{
		return _data.size() ;
	}

	// Get min and max input space values
	vec min() const 
	{
		vec res(4);
		res[0] = 0.0 ;
		res[1] = 0.0 ;
		res[2] = 0.0 ;
		res[3] = 0.0 ;
		return res ;
	}
	vec max() const
	{
		vec res(4);
		res[0] = M_PI / 2 ;
		res[1] = 0.0 ;
		res[2] = M_PI / 2 ;
		res[3] = M_PI ;
		return res ;
	}
};

ALTA_DLL_EXPORT data* provide_data()
{
    return new ASTM();
}

