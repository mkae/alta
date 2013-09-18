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
    return new spherical_gaussian_function();
}

// Overload the function operator
vec spherical_gaussian_function::operator()(const vec& x) const 
{
	return value(x);
}
vec spherical_gaussian_function::value(const vec& x) const 
{
	vec res(dimY());
	double dot;
	params::convert(&x[0], params::CARTESIAN, params::COS_TH, &dot);

	for(int i=0; i<dimY(); ++i)
	{
		res[i] = _ks[i] * std::exp(_n[i] * (dot-1));
	}

	return res;
}

// Reset the output dimension
void spherical_gaussian_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _n.resize(_nY) ;
    _ks.resize(_nY) ;
}

//! Number of parameters to this non-linear function
int spherical_gaussian_function::nbParameters() const 
{
    return 2*dimY();
}

//! Get the vector of parameters for the function
vec spherical_gaussian_function::parameters() const 
{
    vec res(2*dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i*2 + 0] = _ks[i];
        res[i*2 + 1] = _n[i];
    }

    return res;
}

//! \brief get the min values for the parameters
vec spherical_gaussian_function::getParametersMin() const
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
void spherical_gaussian_function::setParameters(const vec& p) 
{
    for(int i=0; i<dimY(); ++i)
    {
        _ks[i] = p[i*2 + 0];
        _n[i]  = p[i*2 + 1];
    }
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec spherical_gaussian_function::parametersJacobian(const vec& x) const 
{
	double dot;
	params::convert(&x[0], params::CARTESIAN, params::COS_TH, &dot);

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j)
			 {

				 // df / dk_s
				 jac[i*nbParameters() + j*2+0] = exp(_n[j] * (dot - 1));

				 // df / dN
				 jac[i*nbParameters() + j*2+1] = _ks[j] * (dot-1) * exp(_n[j]*(dot - 1));
			 }
			 else
			 {
				 jac[i*nbParameters() + j*2+0] = 0.0;
				 jac[i*nbParameters() + j*2+1] = 0.0;
			 }
		 }

    return jac;
}
		
void spherical_gaussian_function::bootstrap(const data* d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = 1.0;
		_n[i]  = 1.0;
	}
}

//! Load function specific files
void spherical_gaussian_function::load(std::istream& in)
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
   if(token.compare("nonlinear_function_spherical_gaussian") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
	}

	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _ks[i];
		in >> token >> _n[i];
	}
}


void spherical_gaussian_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
		out << "#FUNC nonlinear_function_spherical_gaussian" << std::endl ;

		 for(int i=0; i<_nY; ++i)
		 {
			 out << "Ks " << _ks[i] << std::endl;
			 out << "N  " <<  _n[i] << std::endl;
		 }

		 out << std::endl;
	 }
	 else
	 {
		 out << "spherical_gaussian(L, V, N, X, Y, vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _ks[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _n[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "))";
	 }

}

void spherical_gaussian_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 spherical_gaussian(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 Nl)" << std::endl;
        out << "{" << std::endl;
		  out << "\tvec3 H = normalize(L + V);" << std::endl;
        out << "\tvec3 ext_dot = vec3(dot(H,N));" << std::endl;
        out << "\treturn ks * exp(Nl * (ext_dot - vec3(1)));" << std::endl;
        out << "}" << std::endl;
    }
}
