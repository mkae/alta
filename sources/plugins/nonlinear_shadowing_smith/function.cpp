/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <core/common.h>

using namespace alta;

ALTA_DLL_EXPORT function* provide_function()
{
	return new smith();
}

//! Load function specific files
bool smith::load(std::istream& in)
{
//    fresnel::load(in);

    // Parse line until the next comment
    while(in.peek() != '#')
    {
        char line[256];
        in.getline(line, 256);

		  // If we cross the end of the file, or the badbit is
		  // set, the file cannot be loaded
		  if(!in.good())
			  return false;
    }

    // Checking for the comment line #FUNC nonlinear_fresnel_smith
    std::string token;
    in >> token;
    if(token != "#FUNC")
    {
        std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
        return false;
    }

    in >> token;
    if(token != "nonlinear_fresnel_smith")
    {
        std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
        return false;
    }

    // R [double]
    for(int i=0; i<dimY(); ++i)
    {
        in >> token >> w[i];
    }
    return true;
}

void smith::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
        out << "#FUNC nonlinear_fresnel_smith" << std::endl ;
        for(int i=0; i<dimY(); ++i)
        {
            out << "R " << w[i] << std::endl;
        }
        out << std::endl;
    }
    else
    {
        out << "shadowing_smith(L, V, N, X, Y, vec3(";
        for(int i=0; i<dimY(); ++i)
        {
            out << w[i];
            if(i < _nY-1) { out << ", "; }
        }
        out << "))";
    }
}

void smith::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        std::cerr << " FINISH IMPLEMENTATION AT " << __FILE__ << " " << __LINE__ << std::endl;
        assert(0);

        // Generate code for erf function
        // This code is the same as the one in common.hpp 
        out << "/** The Error function in GLSL **/" << std::endl;
        out << "float erf( float a_x ) " << std::endl
            << "{" << std::endl
            << "\t const float a1 = 0.254829592;" << std::endl
            << "\t const float a2 = -0.284496736;" << std::endl
            << "\t const float a3 = 1.421413741;" << std::endl
            << "\t const float a4 = -1.453152027;" << std::endl
            << "\t const float a5 = 1.061405429;" << std::endl
            << "\t const float  p = 0.3275911;" << std::endl
            << "\t // Save the sign of a_x " << std::endl
            << "\t int sign = 1;" << std::endl
            << "\t if( a_x < 0 ) { sign = -1; }" << std::endl
            << "\t const float x = abs( a_x); " << std::endl
            << "\t const float t = 1.0 / (1.0 + p*x); " << std::endl
            << "\t const float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);" << std::endl
            << "\t return sign * y" << std::endl
            
            << "}" << std::endl;

        out << std::endl;

        //Now the shadowing term of smith 
        out << "vec3 shadowing_smith(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 R)" << std::endl;
        out << "{" << std::endl;
        out << "\t const float SQRT_ONE_OVER_TWO_PI = 0.398942280401439;" << std::endl;
        out << "\t const float SQRT_2_OVER_2 = 0.70710678118655;" << std::endl;
        out << "}" << std::endl;
        out << std::endl;        
    
    }
}

vec smith::value(const vec& x) const
{
	vec res(dimY());
    //RP: The Step functions are missing ???? !!!
	if(x[2] == 1.0)
	{
		for(int i=0; i<dimY(); ++i)
		{
			res[i] = 1.0;
		}
	}
	else
	{
        // QUESTION from RP. 
        // This seems to be the monodirectional shadowing term from Smith
        // But in Walter2007 there is an explicit mention to the product
        // of two monodirectional shadowing term G(i,o,h) = G(i,h) * G(v,h)
		double const mu = x[5] / sqrt(x[3]*x[3] + x[4]*x[4]);
    
        double const SQRT_ONE_OVER_TWO_PI = 0.398942280401439 ;
        double const SQRT_2_OVER_2        = 0.70710678118655  ;
		for(int i=0; i<dimY(); ++i)
		{
			double const r = mu / w[i];
			//double const A = sqrt(0.5 / M_PI) * exp(- 0.5 * (r*r)) / r - 0.5 * (1.0 - erf(sqrt(0.5) * r));
            double const A =  SQRT_ONE_OVER_TWO_PI * exp(- 0.5 * (r*r)) / r - 0.5 * (1.0 - erf( SQRT_2_OVER_2 * r));
            
			res[i] = 1.0 / (1.0 + A);
		}

        //RP: In my opinion this code should be add
        // vec res2( dimY() );
        // double const mu2 = x[2] / sqrt(x[0]*x[0] + x[1]*x[1]);

        // for(int i=0; i<dimY(); ++i)
        // {
        //     double const r = mu2 / w[i];
        //     double const A = sqrt(0.5 / M_PI) * exp(- 0.5 * (r*r)) / r - 0.5 * (1.0 - erf(sqrt(0.5) * r));
        //     res2[i] = 1.0 / (1.0 + A);
        // }

        // return res.cwiseProduct( res2 );

	}
	return res;
}

//! \brief Number of parameters to this non-linear function
int smith::nbParameters() const
{
    return dimY();
}

vec smith::getParametersMin() const
{
    vec m(dimY());
    for(int i=0; i<dimY(); ++i) { m[i] = 0.0; }
    return m;
}

//! \brief Get the vector of parameters for the function
vec smith::parameters() const
{
    vec p(dimY());
    for(int i=0; i<dimY(); ++i) { p[i] = w[i]; }
    return p;
}

//! \brief Update the vector of parameters for the function
void smith::setParameters(const vec& p)
{
    for(int i=0; i<dimY(); ++i) { w[i] = p[i]; }
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters.
vec smith::parametersJacobian(const vec& x) const
{
    const int nY = dimY();
    vec jac(nY*nY);

	 vec v = value(x);

    for(int i=0; i<nY; ++i)
        for(int j=0; j<nY; ++j)
        {
			  if(i == j && x[2] < 1.0) 
			  {
				  const double mu = x[5] / sqrt(x[3]*x[3] + x[4]*x[4]);
				  const double r = mu / w[i];

				  jac[j*dimY() + i] = (-1.0 / pow(1.0 + v[i], 2)) * exp(-0.5 * r*r) * (r*r*r/w[i] - 2.0*r/w[i]) / sqrt(2.0 * M_PI) ;
			  }
			  else 
			  {
				  jac[j*dimY() + i] = 0.0;
			  }
		  }

    return jac;
}


void smith::bootstrap(const ptr<data> d, const arguments& args)
{
    for(int i=0; i<dimY(); ++i) { w[i] = 1.0; }
}
