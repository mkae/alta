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
    return new isotropic_lafortune_function();
}

// Overload the function operator
vec isotropic_lafortune_function::operator()(const vec& x) const 
{
	return value(x);
}
vec isotropic_lafortune_function::value(const vec& x) const 
{
	// Get BRDF value
	vec res(dimY());

	double dx, dy, dz;
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];

	// For each color channel
	for(int i=0; i<dimY(); ++i)
	{
		res[i] = _kd[i];

		// For each lobe
		for(int n=0; n<_n; ++n)
		{
			double Cx, Cz, N;
			getCurrentLobe(n, i, Cx, Cz, N);
			double d = Cx*(dx + dy) + Cz*dz;

			if(d > 0.0)
			{
				res[i] += pow(d, N);
			}
		}

#ifndef DEBUG
		if(isnan(res[i]) || res[i] == std::numeric_limits<double>::infinity())
		{
			std::cout << "<<ERROR>> invalid value for input: " << x << std::endl;
		}
#endif	
	}

	return res;
}
        
vec isotropic_lafortune_function::value(const vec& x, const vec& p) const
{
	// Test input parameters for correct size
	assert(p.size() == nbParameters());

	// Get BRDF value
	vec res(dimY());

	double dx, dy, dz;
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];

	// For each color channel
	for(int i=0; i<dimY(); ++i)
	{
		// Start with the diffuse term
		res[i] = _kd[i];

		// For each lobe
		for(int n=0; n<_n; ++n)
		{
			double Cx, Cz, N;
			Cx = p[(n*dimY() + i)*3 + 0];
			Cz = p[(n*dimY() + i)*3 + 1];
			N  = p[(n*dimY() + i)*3 + 2];

			const double d = Cx*(dx + dy) + Cz*dz;
			if(d > 0.0)
			{
				res[i] += pow(d, N);
			}
		}

	}

	return res;
}

// Set the number of lobes of the Lafortune BRDF
void isotropic_lafortune_function::setNbLobes(int N)
{
    _n = N;

    // Update the length of the vectors
    _C.resize(_n*_nY*2) ;
    _N.resize(_n*_nY) ;
}

// Reset the output dimension
void isotropic_lafortune_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _C.resize(_n*_nY*2) ;
    _N.resize(_n*_nY) ;
    _kd.resize(_nY);

    for(int i=0; i<nY; ++i)
        _kd[i] = 0.0;
}

//! Number of parameters to this non-linear function
int isotropic_lafortune_function::nbParameters() const 
{
    return (3*_n)*dimY();
}

//! Get the vector of parameters for the function
vec isotropic_lafortune_function::parameters() const 
{
    vec res((3*_n)*dimY());
    for(int n=0; n<_n; ++n)
	    for(int i=0; i<dimY(); ++i)
		 {
			  res[(n*dimY() + i)*3 + 0] = _C[(n*dimY() + i)*2 + 0];
			  res[(n*dimY() + i)*3 + 1] = _C[(n*dimY() + i)*2 + 1];
			  res[(n*dimY() + i)*3 + 2] = _N[n*dimY()  + i];
		 }

    return res;
}

//! \brief get the min values for the parameters
vec isotropic_lafortune_function::getParametersMin() const
{
    vec res((3*_n)*dimY());
    for(int n=0; n<_n; ++n)
	    for(int i=0; i<dimY(); ++i)
		 {
			  res[(n*dimY() + i)*3 + 0] = -std::numeric_limits<double>::max();
			  res[(n*dimY() + i)*3 + 1] = -std::numeric_limits<double>::max();
			  res[(n*dimY() + i)*3 + 2] = 0.0;
		 }

    return res;
}

//! Update the vector of parameters for the function
void isotropic_lafortune_function::setParameters(const vec& p) 
{
	// Safety check the number of parameters
	assert(p.size() == nbParameters());

	for(int n=0; n<_n; ++n)
		for(int i=0; i<dimY(); ++i)
		{
			_C[(n*dimY() + i)*2 + 0] = p[(n*dimY() + i)*3 + 0];
			_C[(n*dimY() + i)*2 + 1] = p[(n*dimY() + i)*3 + 1];
			_N[n*dimY()  + i]        = p[(n*dimY() + i)*3 + 2];
		}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec isotropic_lafortune_function::parametersJacobian(const vec& x) const 
{
	double dx, dy, dz;
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int n=0; n<_n; ++n)
			 for(int j=0; j<dimY(); ++j)
			 {
				 // index of the current monochromatic lobe
				 int index = i*nbParameters() + 3*(n*dimY() + j);
				 
				 double Cx, Cz, N;
				 getCurrentLobe(n, j, Cx, Cz, N);
				 
				 double d  = Cx*(dx + dy) + Cz*dz;

				 if(i == j && d > 0.0)
				 {
					 // df / dCx
					 jac[index+0] = (dx + dy) * N * std::pow(d, N-1.0);
					 
					 // df / dCz
					 jac[index+1] = dz * N * std::pow(d, N-1.0);

					 // df / dN
					 if(d <= 0.0)
						 jac[index+2] = 0.0;
					 else
						 jac[index+2] = std::log(d) * std::pow(d, N);
				 }
				 else
				 {
					 jac[index+0] = 0.0;
					 jac[index+1] = 0.0;
					 jac[index+2] = 0.0;
				 }
			 }
	 }

    return jac;
}
		
void isotropic_lafortune_function::bootstrap(const data* d, const arguments& args)
{
    // Check the arguments for the number of lobes
    this->setNbLobes(args.get_int("lobes", 1));

    // Set the diffuse component
	vec x0 = d->get(0);
	for(int i=0; i<d->dimY(); ++i)
		_kd[i] = x0[d->dimX() + i];

	for(int i=1; i<d->size(); ++i)
	{
		vec xi = d->get(i);
		for(int j=0; j<d->dimY(); ++j)
			_kd[j] = std::min(xi[d->dimX() + j], _kd[j]);
	}
	std::cout << "<<INFO>> found diffuse: " << _kd << std::endl;

	// Upon user request, the starting position of the lobe can be either load
	// from a file, a distribution beetwen forward backward and dot directions,
	// etc.
	if(args.is_defined("bootstrap"))
	{
		for(int n=0; n<_n; ++n)
		{
			double Cxy, Cz;
			int mod = n % 3;
			if(mod == 0)
			{
				Cxy = -1;
				Cz  =  1;
			}
			else if(mod ==1)
			{
				Cxy = 1;
				Cz  = 1;
			}
			else
			{
				Cxy = 0;
				Cz  = 1;
			}

			for(int i=0; i<dimY(); ++i)
			{
				_C[(n*dimY() + i)*2 + 0] = Cxy;
				_C[(n*dimY() + i)*2 + 1] = Cz;
				_N[n*dimY()  + i]        = (double)_n;
			}
		}
	}
	// The default behaviour of the isotropic Lafortune distribution is to
	// set the transform to [0,0,1] and vary the exponent.
	else
	{
		for(int n=0; n<_n; ++n)
			for(int i=0; i<dimY(); ++i)
			{
				_C[(n*dimY() + i)*2 + 0] = 0.0;
				_C[(n*dimY() + i)*2 + 1] = 1.0;
				_N[n*dimY()  + i]        = (double)_n;
			}
	}
}

std::ofstream& type_definition(std::ofstream& out, int nY)
{
    if(nY == 1)
        out << "float " ;
    else
        out << "vec" << nY ;

    return out;
}

std::ofstream& type_affectation(std::ofstream& out, const std::string& name, const vec& x, int  nY, int n=0, int s=0, int S=1)
{

    out << name << " = ";

	if(nY != 1)
		out << "vec" << nY << "(";

	for(int i=0; i<nY; ++i)
	{
		if(i != 0) out << ", ";
        out << x[n*nY*S + i*S+s];
	}

	if(nY != 1)
		out << ")";

	out << ";" << std::endl;

	return out;
}

//! Load function specific files
void isotropic_lafortune_function::load(std::istream& in)
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
   if(token.compare("nonlinear_function_lafortune") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
	}

	// Shoudl have the #NB_LOBES [int]
	int nb_lobes;
	in >> token >> nb_lobes;
	setNbLobes(nb_lobes);

	// Parse the lobe
	for(int n=0; n<_n; ++n)
	{
		for(int i=0; i<_nY; ++i)
		{

			in >> token >> _C[(n*_nY + i)*2 + 0];
			in >> token >> _C[(n*_nY + i)*2 + 1];
			in >> token >> _N[i];
		}

	}

	std::cout << "<<INFO>> Cd = " << _C << std::endl;
	std::cout << "<<INFO>> N = " << _N << std::endl;

}


void isotropic_lafortune_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
        out << "#FUNC nonlinear_function_lafortune" << std::endl ;
        out << "#NB_LOBES " << _n << std::endl ;

        for(int n=0; n<_n; ++n)
        {

            for(int i=0; i<_nY; ++i)
            {
                out << "Cxy " << _C[(n*_nY + i)*2 + 0] << std::endl;
                out << "Cz  " << _C[(n*_nY + i)*2 + 1] << std::endl;
                out << "N   " << _N[n*_nY + i] << std::endl;

            }


            out << std::endl;
        }
    }
    else
    {
        for(int n=0; n<_n; ++n)
        {
            out << "lafortune(L, V, N, X, Y, vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _C[(n*_nY + i)*2 + 0];
                if(i < _nY-1) { out << ", "; }
            }

            out << "), vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _C[(n*_nY + i)*2 + 1];
                if(i < _nY-1) { out << ", "; }
            }

            out << "), vec3(";
            for(int i=0; i<_nY; ++i)
            {
                out << _N[n*_nY + i];
                if(i < _nY-1) { out << ", "; }
            }

            // For multiple lobes, add a sum sign
            out << "))";
            if(n < _n-1) { out << " + "; }
        }
    }
}

void isotropic_lafortune_function::save_body(std::ostream& out, const arguments& args) const
{
    out << "vec3 lafortune(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 Cx, vec3 Cz, vec3 Nl)" << std::endl;
    out << "{" << std::endl;
    out << "\tvec3 ext_dot = Cx * (dot(L,X)*dot(V,X) + dot(L,Y)*dot(V,Y)) + Cz * dot(L,N)*dot(V,N);" << std::endl;
    out << "\treturn pow(max(ext_dot, vec3(0,0,0)), Nl);" << std::endl;
    out << "}" << std::endl;
}
