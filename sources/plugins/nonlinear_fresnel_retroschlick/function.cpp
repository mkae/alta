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

vec retro_schlick::fresnelValue(const vec& x) const
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
int retro_schlick::nbFresnelParameters() const 
{
    return dimY();
}

//! \brief Get the vector of parameters for the function
vec retro_schlick::getFresnelParameters() const 
{
    vec p(dimY());
    for(int i=0; i<dimY(); ++i) { p[i] = R[i]; }
	return p;
}

//! \brief Update the vector of parameters for the function
void retro_schlick::setFresnelParameters(const vec& p) 
{
    for(int i=0; i<dimY(); ++i) { R[i] = p[i]; }
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters. 
vec retro_schlick::getFresnelParametersJacobian(const vec& x) const 
{
    const int nY = dimY();
    double xp[3], cart[6];
    params::convert(&x[0], input_parametrization(), params::SCHLICK_VK, xp);
    params::convert(&x[0], input_parametrization(), params::CARTESIAN, cart);

    const double dotRK = -xp[0]*cart[0] -xp[1]*cart[1] + xp[2]*cart[2];

	vec jac(nY);
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


void retro_schlick::fresnelBootstrap(const data* d, const arguments& args)
{
    for(int i=0; i<dimY(); ++i) {R[i] = 0.5; }
}

//! Load function specific files
void retro_schlick::load(std::istream& in)
{
	fresnel::load(in);

	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);
	}

	// Checking for the comment line #FUNC nonlinear_fresnel_retro_schlick
	std::string token;
	in >> token;
	if(token != "#FUNC") { std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; }

	in >> token;
	if(token != "nonlinear_fresnel_retro_schlick") { std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; }

	// R [double]
    for(int i=0; i<dimY(); ++i) {
        in >> token >> R[i];
    }
}

void retro_schlick::save_call(std::ostream& out, const arguments& args) const
{
	bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
       f->save_call(out, args);
    }
    else
	{
		out << "(";	f->save_call(out, args); out << ")";
	}


	if(is_alta)
	{
		out << "#FUNC nonlinear_fresnel_retro_schlick" << std::endl ;
		out << "R " << R << std::endl;
		out << std::endl;
	}
	else
	{
		out << " * retro_schlick_fresnel(L, V, N, X, Y, " << R << ")";
	}
}

void retro_schlick::save_body(std::ostream& out, const arguments& args) const
{
	f->save_body(out, args);
	bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

	if(is_shader)
	{
		out << std::endl;
		out << "vec3 retro_schlick_fresnel(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, float R)" << std::endl;
		out << "{" << std::endl;
        out << "\tvec3 R = 2.0f*dot(V,N)*N - V;" << std::endl;
		out << "\tvec3 K = normalize(L + R);" << std::endl;
		out << "\treturn vec3(R + (1.0f - R) * pow(1.0f - clamp(dot(K,R), 0.0f, 1.0f), 5));" << std::endl;
		out << "}" << std::endl;
        out << std::endl;
    }

}


