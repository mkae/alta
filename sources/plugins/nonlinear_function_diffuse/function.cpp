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
    return new diffuse_function();
}

// Overload the function operator
vec diffuse_function::operator()(const vec& x) const 
{
	return value(x);
}
vec diffuse_function::value(const vec& x) const 
{
    vec res(dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i] = _kd[i];
    }

    return res;
}

//! Load function specific files
bool diffuse_function::load(std::istream &in)
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

    // Checking for the comment line #FUNC nonlinear_function_diffuse
    std::string token;
    in >> token;
    if(token != "#FUNC")
    {
        std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl;
        return false;
    }

    in >> token;
    if(token != "nonlinear_function_diffuse")
    {
        std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
        return false;
    }

    // kd [double]
    for(int i=0; i<dimY(); ++i)
    {
        in >> token >> _kd[i];
    }
    return true;

}

void diffuse_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
        out << "#FUNC nonlinear_function_diffuse" << std::endl ;
        for(int i=0; i<dimY(); ++i)
        {
            out << "kd " << _kd[i] << std::endl;
        }
        out << std::endl;
    }
    else
    {
        out << "vec3(";
        for(int i=0; i<dimY(); ++i)
        {
            out << _kd[i]; if(i < dimY()-1) { out << ", "; }
        }
        out << ")";
    }

}
		
//! Number of parameters to this non-linear function
int diffuse_function::nbParameters() const 
{
#ifdef FIT_DIFFUSE
    return dimY();
#else
    return 0;
#endif
}

//! Get the vector of parameters for the function
vec diffuse_function::parameters() const 
{
#ifdef FIT_DIFFUSE
    vec res(dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i*3 + 0] = _kd[i];
    }
#else
    vec res(1);
#endif
    return res;
}

//! Update the vector of parameters for the function
void diffuse_function::setParameters(const vec& p) 
{
#ifdef FIT_DIFFUSE
    for(int i=0; i<dimY(); ++i)
    {
        _kd[i] = p[i];
    }
#endif
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec diffuse_function::parametersJacobian(const vec& x) const 
{
#ifdef FIT_DIFFUSE
	vec jac(dimY());
	for(int i=0; i<dimY(); ++i)
		for(int j=0; j<dimY(); ++j)
		{
			if(i == j)
			{
				// df / dk_d
				jac[i*dimY() + j] = 1.0;

			}
			else
			{
				jac[i*dimY() + j] = 0.0;
			}
		}
#else
	vec jac(1);
#endif

    return jac;
}


void diffuse_function::bootstrap(const data* d, const arguments& args)
{/*
	// Set the diffuse component
	if(params::is_cosine_weighted(d->output_parametrization()))
	{
		vec cart(6);

		for(int i=0; i<d->dimY(); ++i)
			_kd[i] = std::numeric_limits<double>::max();

		for(int i=1; i<d->size(); ++i)
		{
			vec x = d->get(i);
			params::convert(&x[0], d->input_parametrization(), params::CARTESIAN, &cart[0]);
			double cosine = (cart[2] > 0.0 ? cart[2] : 0.0) * (cart[5] > 0.0 ? cart[5] : 0.0);

			if(cosine > 0.0)
			{
				for(int j=0; j<d->dimY(); ++j)
				{
					_kd[j] = std::min(x[d->dimX() + j] / cosine, _kd[j]);
				}
			}
		}
	}
	else
		*/
	{
		for(int i=0; i<d->dimY(); ++i)
			_kd[i] = std::numeric_limits<double>::max();

		for(int i=0; i<d->size(); ++i)
		{
			vec xi = d->get(i);
			for(int j=0; j<d->dimY(); ++j)
				_kd[j] = std::min(xi[d->dimX() + j], _kd[j]);
		}
	}
	std::cout << "<<INFO>> found diffuse: " << _kd << std::endl;

}

