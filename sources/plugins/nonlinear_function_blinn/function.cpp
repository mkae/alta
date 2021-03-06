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

ALTA_DLL_EXPORT function* provide_function(const alta::parameters& params)
{
    return new blinn_function(params);
}

blinn_function::blinn_function(const alta::parameters& params)
    : nonlinear_function(params.set_input(1, params::COS_TH))
{
    // Update the length of the vectors
    _ks.resize(_parameters.dimY()) ;
    _N.resize(_parameters.dimY()) ;
}


// Overload the function operator
vec blinn_function::operator()(const vec& x) const 
{
	return value(x);
}
vec blinn_function::value(const vec& x) const 
{
    vec res(_parameters.dimY());
    for(int i=0; i<_parameters.dimY(); ++i)
    {
		 // Check if the cosine is below the hoziron
		 if(x[0] > 0.0)
		 {
			 res[i] = _ks[i] * std::pow(x[0], _N[i]);
		 }
		 else
		 {
			 res[i] = 0.0;
		 }
    }

    return res;
}

//! Load function specific files
bool blinn_function::load(std::istream& in)
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

    // Checking for the comment line #FUNC nonlinear_function_blinn
    std::string token;
    in >> token;
    if(token != "#FUNC")
    {
        std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
#ifndef DEBUG
        std::cerr << "<<ERROR>> got \"" << token << "\"" << std::endl;
#endif
        return false;
    }

    in >> token;
    if(token != "nonlinear_function_blinn")
    {
        std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
#ifndef DEBUG
        std::cerr << "<<ERROR>> got \"" << token << "\"" << std::endl;
#endif
        return false;
    }

    // ks [double]
	 // N  [double]
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        in >> token >> _ks[i];
        in >> token >>  _N[i];
    }

#ifdef DEBUG
    std::cout << "<<DEBUG>> load parameters " << parameters() << std::endl;
#endif
    return true;
}

//! Number of parameters to this non-linear function
int blinn_function::nbParameters() const 
{
    return 2*_parameters.dimY();
}

//! Get the vector of parameters for the function
vec blinn_function::parameters() const 
{
    vec res(2*_parameters.dimY());
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        res[i*2 + 0] = _ks[i];
        res[i*2 + 1] = _N[i];
    }

    return res;
}

//! Update the vector of parameters for the function
void blinn_function::setParameters(const vec& p) 
{
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        _ks[i] = p[i*2 + 0];
        _N[i]  = p[i*2 + 1];
    }
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec blinn_function::parametersJacobian(const vec& x) const 
{
    vec jac(_parameters.dimY()*nbParameters());
	 for(int i=0; i<_parameters.dimY(); ++i)
		 for(int j=0; j<_parameters.dimY(); ++j)
		 {
			 if(i == j)
			 {
				 // Test if the configuration is below the horizon
				 if(x[0] <= 0.0)
				 {
					jac[i*nbParameters() + j*2+0] = 0.0;
					jac[i*nbParameters() + j*2+1] = 0.0;
				 }
				 else
				 {
					 // df / dk_s
					 jac[i*nbParameters() + j*2+0] = std::pow(x[0], _N[j]);

					 // df / dN
					 jac[i*nbParameters() + j*2+1] = _ks[j] * log(x[0]) * std::pow(x[0], _N[j]);
				 }
			 }
			 else
			 {
				 jac[i*nbParameters() + j*2+0] = 0.0;
				 jac[i*nbParameters() + j*2+1] = 0.0;
			 }
		 }

    return jac;
}


void blinn_function::bootstrap(const ptr<data> d, const arguments& args)
{
    for(int i=0; i<_parameters.dimY(); ++i)
    {
        _ks[i] = 1.0;
        _N[i]  = 10.0;
    }
}

void blinn_function::save_call(std::ostream& out, 
                               const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
		out << "#FUNC nonlinear_function_blinn" << std::endl ;

		 for(int i=0; i<_parameters.dimY(); ++i)
		 {
			 out << "Ks " << _ks[i] << std::endl;
			 out << "N  " <<  _N[i] << std::endl;
		 }

		 out << std::endl;
	 }
	 else
	 {
		 out << "blinn(L, V, N, X, Y, vec3(";
		 for(int i=0; i<_parameters.dimY(); ++i)
		 {
			 out << _ks[i];
			 if(i < _parameters.dimY()-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_parameters.dimY(); ++i)
		 {
			 out << _N[i];
			 if(i < _parameters.dimY()-1) { out << ", "; }
		 }

		 out << "))";
	 }

}

void blinn_function::save_body(std::ostream& out, 
                               const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 blinn(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 Nl)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3 H = normalize(L + V);" << std::endl;
        out << "\tvec3 ext_dot = vec3(dot(H,N));" << std::endl;
        out << "\treturn ks * pow(max(ext_dot, vec3(0,0,0)), Nl);" << std::endl;
        out << "}" << std::endl;
        out << std::endl;
    }

}


