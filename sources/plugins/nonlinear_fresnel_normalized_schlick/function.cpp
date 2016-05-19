/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

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

#include <core/common.h>

using namespace alta;

ALTA_DLL_EXPORT function* provide_function()
{
	return new schlick();
}

schlick::schlick():
    // XXX: Partial parametrization.
    nonlinear_function(alta::parameters(6, 0,
                                        params::CARTESIAN, params::UNKNOWN_OUTPUT))
{
}

//! Load function specific files
bool schlick::load(std::istream& in)
{
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

	// Checking for the comment line #FUNC nonlinear_fresnel_schlick
	std::string token;
	in >> token;
    if(token != "#FUNC")
    {
        std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
        return false;
    }

	in >> token;
    if(token != "nonlinear_fresnel_normalized_schlick")
    {
        std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
        return false;
    }

	// R [double]
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        in >> token >> R[i];
    }
    return true;
}

void schlick::save_call(std::ostream& out, const arguments& args) const
{
	bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

	if(is_alta)
	{
        out << "#FUNC nonlinear_fresnel_normalized_schlick" << std::endl ;
        for(int i=0; i<_parameters.dimY(); ++i)
        {
            out << "R " << R[i] << std::endl;
        }
		out << std::endl;
	}
	else
	{
        out << "normalized_schlick_fresnel(L, V, N, X, Y, vec3(";
        for(int i=0; i<_parameters.dimY(); ++i)
        {
            out << R[i];
            if(i < _parameters.dimY()-1) { out << ", "; }
        }
        out << "))";
	}
}

void schlick::save_body(std::ostream& out, const arguments& args) const
{
	bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

	if(is_shader)
	{
		out << std::endl;
        out << "vec3 normalized_schlick_fresnel(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 R)" << std::endl;
		out << "{" << std::endl;
        out << "\tvec3 H = normalize(L + V);" << std::endl;
        out << "\treturn vec3(1.0f) + (vec3(1.0f)/R - vec3(1.0f)) * pow(1.0f - clamp(dot(H,V), 0.0f, 1.0f), 5);" << std::endl;
		out << "}" << std::endl;
        out << std::endl;
	}

}

vec schlick::value(const vec& x) const
{
    vec xp(3), cart(6);
    params::convert(&x[0], _parameters.input_parametrization(),
                    params::RUSIN_VH, &xp[0]);
    params::convert(&x[0], _parameters.input_parametrization(),
                    params::CARTESIAN, &cart[0]);
	 //std::cout << xp << " ===== " << cart << std::endl;


    const double dotVH = xp[0]*cart[0] + xp[1]*cart[1] + xp[2]*cart[2];

    vec res(_parameters.dimY());
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        res[i] = 1.0 + (1.0/R[i] - 1.0) * pow(1.0 - clamp(dotVH, 0.0, 1.0), 5.0);
    }
/*
	 std::cout << " # " << R << ", " << dotVH << std::endl;
	 std::cout << " ## " << res << std::endl;
    
	 */
	 return res;
}

//! \brief Number of parameters to this non-linear function
int schlick::nbParameters() const 
{
    return _parameters.dimY();
}

vec schlick::getParametersMin() const
{
    vec m(nbParameters());
    for(int i=0; i<nbParameters(); ++i) { m[i] = 1.0E-10; }
    return m;
}

vec schlick::getParametersMax() const
{
    vec m(nbParameters());
    for(int i=0; i<nbParameters(); ++i) { m[i] = 1.0; }
    return m;
}

//! \brief Get the vector of parameters for the function
vec schlick::parameters() const 
{
    vec p(_parameters.dimY());
    for(int i=0; i<_parameters.dimY(); ++i) { p[i] = R[i]; }
	return p;
}

//! \brief Update the vector of parameters for the function
void schlick::setParameters(const vec& p) 
{
    for(int i=0; i<_parameters.dimY(); ++i) { R[i] = p[i]; }
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters. 
vec schlick::parametersJacobian(const vec& x) const 
{
	const int nY = _parameters.dimY();
    double xp[3], cart[6];
    params::convert(&x[0], _parameters.input_parametrization(),
                    params::RUSIN_VH, xp);
    params::convert(&x[0], _parameters.input_parametrization(),
                    params::CARTESIAN, cart);

    const double dotVH = xp[0]*cart[0] + xp[1]*cart[1] + xp[2]*cart[2];

    vec jac(nY*nY);
	for(int i=0; i<nY; ++i)
    {
        for(int j=0; j<nY; ++j)
        {
            if(i == j) {
                jac[j*_parameters.dimY() + i] = - (1.0/R[i]) * pow(1.0 - clamp(dotVH, 0.0, 1.0), 5.0);
            } else {
                jac[j*_parameters.dimY() + i] = 0.0;
            }
        }
    }

	return jac;
}


void schlick::bootstrap(const ptr<data> d, const arguments& args)
{
    for(int i=0; i<_parameters.dimY(); ++i) 
    { 
        R[i] = 0.5; 
    }
    
    std::cout << "<<INFO>> Normalized Fresnel with Schlick Approx. Fit will start from " << R << std::endl;
}
