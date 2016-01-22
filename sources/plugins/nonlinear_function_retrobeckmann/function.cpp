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
    return new beckmann_function();
}

// Overload the function operator
vec beckmann_function::operator()(const vec& x) const 
{
	return value(x);
}
vec beckmann_function::value(const vec& x) const 
{
	double dot;
	if(_use_back_param)
	{
		double h[3];
		params::convert(&x[0], params::CARTESIAN, params::SCHLICK_VK, &h[0]);
		dot = h[2];
	}
	else
	{
		dot = std::max(x[0]*x[3] +  x[2]*x[4] +  x[2]*x[5], 0.0);
	}

	vec res = vec::Zero(dimY());

	for(int i=0; i<dimY(); ++i)
	{
		const double a    = _a[i];
		const double a2   = a*a;
		const double dh2  = dot*dot;
		const double expo = exp((dh2 - 1.0) / (a2 * dh2));

		if(dot > 0.0)
		{
			res[i] = _ks[i] / (4.0 /* x[2]*x[5] */* M_PI * a2 * dh2*dh2) * expo;
		}
		else
		{
			res[i] = 0.0; 
		}
	}
	return res;
}

// Reset the output dimension
void beckmann_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _a.resize(_nY) ;
    _ks.resize(_nY) ;
}

//! Number of parameters to this non-linear function
int beckmann_function::nbParameters() const 
{
	return 2*dimY();
}

//! Get the vector of parameters for the function
vec beckmann_function::parameters() const 
{
	vec res(2*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*2 + 0] = _ks[i];
		res[i*2 + 1] = _a[i];
	}
	return res;
}

//! \brief get the min values for the parameters
vec beckmann_function::getParametersMin() const
{
	vec res(2*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*2 + 0] = 0.0;
		res[i*2 + 1] = 0.0;
	}

	return res;
}

//! Update the vector of parameters for the function
void beckmann_function::setParameters(const vec& p) 
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = p[i*2 + 0];
		_a[i]  = p[i*2 + 1];
	}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters
//! \todo finish. 
vec beckmann_function::parametersJacobian(const vec& x) const 
{
	double dot;
	if(_use_back_param)
	{
		double h[3];
		params::convert(&x[0], params::CARTESIAN, params::SCHLICK_VK, &h[0]);
		dot = h[2];
	}
	else
	{
		dot = std::max(x[0]*x[3] +  x[2]*x[4] +  x[2]*x[5], 0.0);
	}

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j && dot>0.0/* && x[2]*x[5]>0.0*/)
			 {
				 const double a    = _a[i];
				 const double a2   = a*a;
				 const double dh2  = dot*dot;
				 const double expo = exp((dh2 - 1.0) / (a2 * dh2));
				 const double fac  = (4.0 /* x[2]*x[5] */* M_PI * a2 * dh2*dh2);

				 // df / dk_s
				 jac[i*nbParameters() + j*2+0] = expo / fac;

				 // df / da_x
				 jac[i*nbParameters() + j*2+1] = - _ks[i] * (expo/(4.0/*x[2]*x[5]*/)) * ((2* a * dot)/(M_PI*a2*a2*dh2)) * (1 + (dh2 - 1.0)*dot/(a2*dh2*dot));
			 }
			 else
			 {
				 jac[i*nbParameters() + j*2+0] = 0.0;
				 jac[i*nbParameters() + j*2+1] = 0.0;
			 }
		 }
	 }

    return jac;
}
		
void beckmann_function::bootstrap(const ptr<data> d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = 1.0;
		_a[i]  = 1.0;
	}

	if(args.is_defined("retro"))
	{
		_use_back_param = false;
	}
	else
	{
		_use_back_param = true;
	}
}

//! Load function specific files
bool beckmann_function::load(std::istream& in)
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
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
		return false;
	}

	in >> token;
   if(token.compare("nonlinear_function_retrobeckmann") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
		return false;
	}
	
	in >> token;
	if(token.compare("#TYPE") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #TYPE is not the next line defined." << std::endl; 
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
		return false;
	}

	in >> token;
	if(token.compare("RETRO") != 0) 
	{ 
		_use_back_param = false;
	}
	else
	{
		_use_back_param = true;
	}


	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _ks[i];
		in >> token >> _a[i];
	}
    return true;
}


void beckmann_function::save_call(std::ostream& out, const arguments& args) const
{
	bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

	if(is_alta)
	{
		out << "#FUNC nonlinear_function_retrobeckmann" << std::endl ;
		out << "#TYPE ";
		if(_use_back_param) { out << "BACK" << std::endl; }
		else { out << "RETRO" << std::endl; }

			for(int i=0; i<_nY; ++i)
			{
				out << "Ks " << _ks[i] << std::endl;
				out << "a  " << _a[i]  << std::endl;
			}

		out << std::endl;
	}
	else
	{
		out << "retrobeckmann(L, V, N, X, Y, vec3(";
		for(int i=0; i<_nY; ++i)
		{
			out << _ks[i];
			if(i < _nY-1) { out << ", "; }
		}

		out << "), vec3(";
		for(int i=0; i<_nY; ++i)
		{
			out << _a[i];
			if(i < _nY-1) { out << ", "; }
		}
		out << "))";
	}
}

void beckmann_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 retrobeckmann(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 a)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3  R   = 2*dot(V,N)*N - V;" << std::endl;
        out << "\tvec3  B   = normalize(L + R);" << std::endl;
        out << "\tfloat bn  = dot(B,N);" << std::endl;
		  out << "\tfloat ln  = dot(L,N);" << std::endl;
		  out << "\tfloat vn  = dot(V,N);" << std::endl;
		  out << "\t" << std::endl;
        out << "\treturn ks / (4 * " << M_PI << " * a*a * ln*vn) * exp((bn*bn - 1.0) / (a*a*bn*bn));" << std::endl;
        out << "}" << std::endl;
    }
}
