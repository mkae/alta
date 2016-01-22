/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

using namespace alta;

ALTA_DLL_EXPORT function* provide_function()
{
    return new abc_function();
}

// Overload the function operator
vec abc_function::operator()(const vec& x) const 
{
	return value(x);
}
vec abc_function::value(const vec& x) const 
{

	// Compute the Shadow term to init res
	vec res = vec::Zero(dimY());
	const double hn = 1.0 - x[0];

	for(int i=0; i<dimY(); ++i)
	{
		if(hn >= 0.0 && hn <= 1.0)
		{
			res[i] = _a[i] / pow(1.0 + _b[i]*hn, _c[i]);
		}
		else
		{
			res[i] = 0.0;
		}
	}
	return res;
}

// Reset the output dimension
void abc_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _a.resize(_nY) ;
    _b.resize(_nY) ;
    _c.resize(_nY) ;
}

//! Number of parameters to this non-linear function
int abc_function::nbParameters() const 
{
	return 3*dimY();
}

//! Get the vector of parameters for the function
vec abc_function::parameters() const 
{
	vec res(nbParameters());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*3 + 0] = _a[i];
		res[i*3 + 1] = _b[i];
		res[i*3 + 2] = _c[i];
	}
	return res;
}

//! \brief get the min values for the parameters
vec abc_function::getParametersMin() const
{
	vec res(nbParameters());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*3 + 0] = 0.0;
		res[i*3 + 1] = 0.0;
		res[i*3 + 2] = 0.0;
	}

	return res;
}

//! Update the vector of parameters for the function
void abc_function::setParameters(const vec& p) 
{
	for(int i=0; i<dimY(); ++i)
	{
		_a[i]  = p[i*3 + 0];
		_b[i]  = p[i*3 + 1];
		_c[i]  = p[i*3 + 2];
	}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters
//! \todo finish. 
vec abc_function::parametersJacobian(const vec& x) const 
{
    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j)
			 {
				 const double hn = 1.0 - x[0];
				 const double f  = 1.0 + _b[i]*hn;
				 const double denom = pow(f, _c[i]);
				 const double fact  = 1.0 / denom;
				 const double fact2 = fact*fact;

				 // df / da
				 jac[i*nbParameters() + j*3+0] = fact;

				 // df / db
				 jac[i*nbParameters() + j*3+1] = - _a[i] * hn *_c[i] * pow(f, _c[i]-1.0) * fact2;

				 // df / dc
				 if(f > 0.0)
				 {
					 jac[i*nbParameters() + j*3+2] = - _a[i] * log(f) * fact;
				 }
				 else
				 {
					jac[i*nbParameters() + j*3+2] = 0.0;
				 }
			 }
			 else
			 {
				 jac[i*nbParameters() + j*3+0] = 0.0;
				 jac[i*nbParameters() + j*3+1] = 0.0;
				 jac[i*nbParameters() + j*3+2] = 0.0;
			 }
		 }
	 }

    return jac;
}
		
void abc_function::bootstrap(const ptr<data> d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i)
	{
		_a[i]  = 1.0;
		_b[i]  = 1.0;
		_c[i]  = 1.0;
	}
	
	if(args.is_defined("param"))
	{
		setParametrization(params::parse_input(args["param"]));
	}
	else
	{
		setParametrization(params::COS_TH);
	}
	if(params::dimension(input_parametrization()) != 1)
	{
		std::cerr << "<<ERROR>> the parametrization specifed in the file for the ABC model is incorrect" << std::endl;
	}
}

//! Load function specific files
bool abc_function::load(std::istream& in)
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

    // Checking for the comment line #FUNC nonlinear_function_lafortune
	std::string token;
	in >> token;
	if(token.compare("#FUNC") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
        return false;
	}

	in >> token;
   if(token.compare("nonlinear_function_abc") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
        return false;
	}
	
	in >> token;
	if(token.compare("#PARAM") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #PARAM is not the next line defined." << std::endl; 
        return false;
	}
	in >> token;
	setParametrization(params::parse_input(token));
	if(params::dimension(input_parametrization()) != 1)
	{
		std::cerr << "<<ERROR>> the parametrization specifed in the file for the ABC model is incorrect" << std::endl;
	}

	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _a[i];
		in >> token >> _b[i];
		in >> token >> _c[i];
	}
    return true;
}


void abc_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
		out << "#FUNC nonlinear_function_abc" << std::endl ;
		out << "#PARAM " << params::get_name(input_parametrization()) << std::endl;

		 for(int i=0; i<_nY; ++i)
		 {
			 out << "a  " << _a[i]  << std::endl;
			 out << "b  " << _b[i]  << std::endl;
			 out << "c  " << _c[i]  << std::endl;
		 }
	
		 out << std::endl;
	 }
	 else
	 {
		 out << "abc(L, V, N, X, Y, vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _a[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _b[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _c[i];
			 if(i < _nY-1) { out << ", "; }
		 }
		 out << "))";
	 }

}

void abc_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 abc(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 a, vec3 b, vec3 c)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3  H   = normalize(L + V);" << std::endl;
        out << "\tfloat hn  = dot(H,N);" << std::endl;
		  out << "\t" << std::endl;
        out << "\treturn a / pow(vec3(1.0f) + b*hn, c);" << std::endl;
        out << "}" << std::endl;
    }
}
