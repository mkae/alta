#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT function* provide_function()
{
	return new retro_schlick();
}

retro_schlick::retro_schlick()
{
	setParametrization(params::CARTESIAN);
	setDimX(6);
}

vec retro_schlick::value(const vec& x) const
{
	double xp[3], cart[6];
	params::convert(&x[0], input_parametrization(), params::SCHLICK_VK, xp);
	params::convert(&x[0], input_parametrization(), params::CARTESIAN, cart);

	const double dotRK = -xp[0]*cart[0] -xp[1]*cart[1] + xp[2]*cart[2];

	vec res(_nY);
	for(int i=0; i<_nY; ++i)
	{
		res[i] = R[i] + (1.0 - R[i]) * pow(1.0 - clamp(dotRK, 0.0, 1.0), 5.0);
	}

	return res;
}

//! \brief Number of parameters to this non-linear function
int retro_schlick::nbParameters() const 
{
	return dimY();
}

//! \brief Get the vector of parameters for the function
vec retro_schlick::parameters() const 
{
	vec p(dimY());
	for(int i=0; i<dimY(); ++i) { p[i] = R[i]; }
	return p;
}

//! \brief Update the vector of parameters for the function
void retro_schlick::setParameters(const vec& p) 
{
	for(int i=0; i<dimY(); ++i) { R[i] = p[i]; }
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters. 
vec retro_schlick::parametersJacobian(const vec& x) const 
{
	const int nY = dimY();
	double xp[3], cart[6];
	params::convert(&x[0], input_parametrization(), params::SCHLICK_VK, xp);
	params::convert(&x[0], input_parametrization(), params::CARTESIAN, cart);

	const double dotRK = -xp[0]*cart[0] -xp[1]*cart[1] + xp[2]*cart[2];

	vec jac(nbParameters()*nY);
	for(int i=0; i<nY; ++i)
		for(int j=0; j<nY; ++j)
		{
			if(i == j)
			{
				jac[i*nY + j] = 1.0 - pow(1.0 - clamp(dotRK, 0.0, 1.0), 5.0);
			}
			else
			{
				jac[i*nY + j] = 0.0;
			}
		}

	return jac;
}

vec retro_schlick::getParametersMin() const
{
	vec m = vec::Zero(dimY());
	return m;
}

vec retro_schlick::getParametersMax() const
{
	vec M(dimY());
	for(int i=0; i<dimY(); ++i) { M[i] = 1.0; }
	return M;
}

void retro_schlick::bootstrap(const data* d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i) {R[i] = 0.5; }
}

//! Load function specific files
bool retro_schlick::load(std::istream& in)
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

	// Checking for the comment line #FUNC nonlinear_fresnel_retro_schlick
	std::string token;
	in >> token;
	if(token != "#FUNC")
	{
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
		return false;
	}

	in >> token;
	if(token != "nonlinear_fresnel_retro_schlick")
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
		return false;
	}

	// R [double]
	for(int i=0; i<dimY(); ++i) {
		in >> token >> R[i];
	}
	return true;
}

void retro_schlick::save_call(std::ostream& out, const arguments& args) const
{
	bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

	if(is_alta)
	{
		out << "#FUNC nonlinear_fresnel_retro_schlick" << std::endl ;
		for(int i=0; i<dimY(); ++i)
		{
			out << "R " << R[i] << std::endl;
		}
		out << std::endl;
	}
	else
	{
		out << "retro_schlick_fresnel(L, V, N, X, Y, vec3(";
		for(int i=0; i<dimY(); ++i)
		{
			out << R[i];
			if(i < _nY-1) { out << ", "; }
		}
		out << "))";
	}
}

void retro_schlick::save_body(std::ostream& out, const arguments& args) const
{
	bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

	if(is_shader)
	{
		out << std::endl;
		out << "vec3 retro_schlick_fresnel(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 R0)" << std::endl;
		out << "{" << std::endl;
		out << "\tvec3 R = 2.0f*dot(V,N)*N - V;" << std::endl;
		out << "\tvec3 K = normalize(L + R);" << std::endl;
		out << "\treturn vec3(R0 + (vec3(1.0f) - R0) * pow(1.0f - clamp(dot(K,R), 0.0f, 1.0f), 5.0f));" << std::endl;
		out << "}" << std::endl;
		out << std::endl;
	}

}


