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
    return new ward_function();
}

// Overload the function operator
vec ward_function::operator()(const vec& x) const 
{
	return value(x);
}
vec ward_function::value(const vec& x) const 
{
	vec res(dimY());

	double h[3];
	params::convert(&x[0], params::CARTESIAN, params::RUSIN_VH, &h[0]);

	for(int i=0; i<dimY(); ++i)
	{
		const double ax = _ax[i];
		const double ay = _ay[i];

		const double hx_ax = h[0]/ax;
		const double hy_ay = h[1]/ay;

		const double exponent = (hx_ax*hx_ax + hy_ay*hy_ay) / (h[2]*h[2]);

		if(x[2]*x[5] > 0.0)
		{
			res[i] = (_ks[i] / (4.0 * M_PI * ax * ay * sqrt(x[2]*x[5]))) * std::exp(- exponent);
		}
		else
		{
			res[i] = 0.0;
		}
	}

	return res;
}

// Reset the output dimension
void ward_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _ax.resize(_nY) ;
    _ay.resize(_nY) ;
    _ks.resize(_nY) ;
}

//! Number of parameters to this non-linear function
int ward_function::nbParameters() const 
{
	return 3*dimY();
}

//! Get the vector of parameters for the function
vec ward_function::parameters() const 
{
	vec res(3*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*3 + 0] = _ks[i];
		res[i*3 + 1] = _ax[i];
        res[i*3 + 2] = _ay[i];
	}
	return res;
}

//! \brief get the min values for the parameters
vec ward_function::getParametersMin() const
{
	vec res(3*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*3 + 0] = 0.0;
		res[i*3 + 1] = 0.0;
		res[i*3 + 2] = 0.0;
	}

	return res;
}

//! Update the vector of parameters for the function
void ward_function::setParameters(const vec& p) 
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = p[i*3 + 0];
		_ax[i] = p[i*3 + 1];
        _ay[i] = p[i*3 + 2];
	}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters
//! \todo finish. 
vec ward_function::parametersJacobian(const vec& x) const 
{
	double h[3];
	params::convert(&x[0], params::CARTESIAN, params::RUSIN_VH, h);

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int j=0; j<dimY(); ++j)
		 {
             if(i == j && x[2]*x[5]>0.0)
			 {
                 const double ax = _ax[i];
                 const double ay = _ay[i];

                 const double hx_ax = h[0]/ax;
                 const double hy_ay = h[1]/ay;
				 const double gauss = exp(-(hx_ax*hx_ax + hy_ay*hy_ay) / (h[2]*h[2]));
                 const double fact  = 1.0 / (4.0*M_PI*ax*ay*sqrt(x[2]*x[5]));

				 // df / dk_s
				 jac[i*nbParameters() + j*3+0] = fact * gauss;

				 // df / da_x
                 jac[i*nbParameters() + j*3+1] = _ks[i] * fact * (1.0/ax) * (((2.0*h[0]*hx_ax) / h[2]) * gauss - 1);

				 // df / da_y
                 jac[i*nbParameters() + j*3+2] = _ks[i] * fact * (1.0/ay) * (((2.0*h[1]*hy_ay) / h[2]) * gauss - 1);
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
		
void ward_function::bootstrap(const ptr<data> d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = 1.0;
		_ax[i] = 1.0;
		_ay[i] = 1.0;
	}

    // Check if the fit should be done on an istropic lobe or not.
    if(args.is_defined("isotropic"))
    {
        isotropic = true;
    }
}

//! Load function specific files
bool ward_function::load(std::istream& in)
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
   if(token.compare("nonlinear_function_ward") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
        return false;
	}

	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _ks[i];
		in >> token >> _ax[i];
		in >> token >> _ay[i];
	}
    return true;
}


void ward_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
		out << "#FUNC nonlinear_function_ward" << std::endl ;

		 for(int i=0; i<_nY; ++i)
		 {
			 out << "Ks " << _ks[i] << std::endl;
			 out << "ax " << _ax[i] << std::endl;
             out << "ay " << _ay[i] << std::endl;
		 }
	
		 out << std::endl;
	 }
	 else
	 {
		 out << "ward(L, V, N, X, Y, vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _ks[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _ax[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _ay[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "))";
	 }

}

void ward_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 ward(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 ax, vec3 ay)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3  H   = normalize(L + V);" << std::endl;
        out << "\tvec3  hax = dot(H,X) / ax;" << std::endl;
        out << "\tvec3  hay = dot(H,Y) / ay;" << std::endl;
        out << "\tfloat hn  = dot(H,N);" << std::endl;
		  out << "\tfloat ln  = dot(L,N);" << std::endl;
		  out << "\tfloat vn  = dot(V,N);" << std::endl;
		  out << "\t" << std::endl;
		  out << "\tif(ln*vn > 0) {" << std::endl;
        out << "\t\treturn (ks / (4 * " << M_PI << " * ax*ay * sqrt(ln*vn))) * exp(-(hax*hax + hay*hay)/(hn*hn));" << std::endl;
		  out << "\t} else {" << std::endl;
		  out << "\t\t return vec3(0);" << std::endl;
		  out << "\t}" << std::endl;
        out << "}" << std::endl;
    }
}
