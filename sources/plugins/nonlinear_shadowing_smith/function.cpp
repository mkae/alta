#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <core/common.h>

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
        // out << std::endl;
        // out << "vec3 shadowing_smith(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 R)" << std::endl;
        // out << "{" << std::endl;
        // out << "\tvec3 H = normalize(L + V);" << std::endl;
        // out << "\treturn R + (vec3(1.0) - R) * pow(1.0f - clamp(dot(H,V), 0.0f, 1.0f), 5);" << std::endl;
        // out << "}" << std::endl;
        // out << std::endl;
        std::cerr << " IMPLEMENT ME AT " << __FILE__ << " " << __LINE__ << std::endl;
        assert(0);
    }

}

vec smith::value(const vec& x) const
{
	vec res(dimY());
	if(x[2] == 1.0)
	{
		for(int i=0; i<dimY(); ++i)
		{
			res[i] = 1.0;
		}
	}
	else
	{
		const double mu = x[5] / sqrt(x[3]*x[3] + x[4]*x[4]);

		for(int i=0; i<dimY(); ++i)
		{
			const double r = mu / w[i];
			const double A = sqrt(0.5 / M_PI) * exp(- 0.5 * (r*r)) / r - 0.5 * (1.0 - erf(sqrt(0.5) * r));
			res[i] = 1.0 / (1.0 + A);
		}
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
