/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "clustering.h"

#include <iostream>
#include <limits>
#include <string>
#include <fstream>

template<class T> void clustering(const T* in_data, int _nY, params::input in_param, params::input out_param, std::vector<vec>& out_data)
{

    // Set the dimensions
    const int in_nX  = params::dimension(in_param);
    const int out_nX = params::dimension(out_param);

#ifdef DEBUG
    std::cout << " in dim = " <<  in_nX << std::endl;
    std::cout << "out dim = " << out_nX << std::endl;
#endif

	 for(int i=0; i<in_data->size(); ++i)
	 {
		 vec p = in_data->get(i);
		 vec e (out_nX + _nY);

		 // Fill the input part of the vector
		 params::convert(&p[0], in_param, out_param, &e[0]);
/*
		 // Search for duplicates only when the reparametrization is compressing
		 // the space.
		 if(out_nX < in_nX)
		 {
			 bool already_exist = false;
			 for(unsigned int j=0; j<out_data.size(); ++j)
			 {
				 vec test_e = out_data[j];

				 double dist = 0.0;
				 for(int k=0; k<out_nX; ++k)
				 {
					 const double temp = test_e[k]-e[k];
					 dist += temp*temp;
				 }

				 if(dist < 1.0E-5)
					 already_exist = true;
			 }

			 if(!already_exist)
			 {
				 // Fill the output part of the vector
				 for(int j=0; j<_nY; ++j)
					 e[out_nX + j] = p[in_nX + j];

#ifdef DEBUG
				 std::cout << " in = " << p << std::endl;
				 std::cout << "out = " << e << std::endl;
#endif
				 out_data.push_back(e);
			 }
		 }
		 else
*/
		 {
			 // Fill the output part of the vector
			 for(int j=0; j<_nY; ++j)
				 e[out_nX + j] = p[in_nX + j];
			 
			 out_data.push_back(e);
		 }
	 }

}
