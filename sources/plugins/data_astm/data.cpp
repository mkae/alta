/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "data.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>


// Load data from a file
void data_astm::load(const std::string& filename) 
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
void data_astm::load(const std::string& filename, const arguments& args)
{
	this->load(filename);
}

// Acces to data
vec data_astm::get(int i) const 
{
	return _data[i];
}
vec data_astm::operator[](int i) const 
{
	return _data[i];
}
vec data_astm::value(vec in, vec out) const
{
	std::cout << "<<ERROR>> not implemented" << std::endl;
	vec res(4);
	return res;
}

// Get data size, e.g. the number of samples to fit
int data_astm::size() const 
{
	return _data.size() ;
}

// Get min and max input space values
vec data_astm::min() const 
{
	vec res(4);
	res[0] = 0.0 ;
	res[1] = 0.0 ;
	res[2] = 0.0 ;
	res[3] = 0.0 ;
	return res ;
}
vec data_astm::max() const
{
	vec res(4);
	res[0] = M_PI / 2 ;
	res[1] = 0.0 ;
	res[2] = M_PI / 2 ;
	res[3] = M_PI ;
	return res ;
}

ALTA_DLL_EXPORT data* provide_data()
{
    return new data_astm();
}

