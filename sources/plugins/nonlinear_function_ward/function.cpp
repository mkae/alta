#include "function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

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
	params::convert(&x[0], params::CARTESIAN, params::RUSIN_VH, h);

	for(int i=0; i<dimY(); ++i)
	{
		const double hx_ax = h[0]/_ax[i];
		const double hy_ay = h[1]/_ay[i];
		const double exponent = (hx_ax*hx_ax + hy_ay*hy_ay) / (h[2]*h[2]);
		res[i] = (_ks[i] / (4.0 * M_PI * _ax[i] * _ay[i] * sqrt(x[2]*x[5]))) * std::exp(- exponent);
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
	double dot = compute_dot(x);

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j)
			 {
/*
				 // df / dk_s
				 jac[i*nbParameters() + j*2+0] = exp(_n[j] * (dot - 1));

				 // df / dN
				 jac[i*nbParameters() + j*2+1] = _ks[j] * (dot-1) * exp(_n[j]*(dot - 1));
*/
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
		
void ward_function::bootstrap(const data* d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = 1.0;
		_ax[i] = 1.0;
		_ay[i] = 1.0;
	}
}

//! Load function specific files
void ward_function::load(std::istream& in)
{
	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);
	}

    // Checking for the comment line #FUNC nonlinear_function_lafortune
	std::string token;
	in >> token;
	if(token.compare("#FUNC") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
	}

	in >> token;
   if(token.compare("nonlinear_function_ward") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
	}

	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _ks[i];
		in >> token >> _ax[i];
		in >> token >> _ay[i];
	}
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
        out << "\treturn (ks / (4 * M_PI * ax*ay * sqrt(dot(L,N)*dot(V,N))) * exp(-(hax*hax + hay*hay)/(hn*hn));" << std::endl;
        out << "}" << std::endl;
    }
}
