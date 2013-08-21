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
				res[i] += pow(d, N);
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
void isotropic_lafortune_function::load(const std::string& filename)
{
	std::ifstream file(filename.c_str()) ;
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
		throw ;
	}

	_nX = 0 ; _nY = 0 ;
	_n = 0;

	double x, y, dy ;
	while(file.peek() == '#')
	{
		std::string line ;
		std::getline(file, line) ;
		std::stringstream linestream(line) ;

		linestream.ignore(1) ;

		std::string comment ;
		linestream >> comment ;

		if(comment == std::string("DIM"))
		{
			linestream >> _nX >> _nY ;
		}
		else if(comment == std::string("NB_LOBES"))
		{
			linestream >> _n ;
		}
	}

	_kd = vec(_nY);
	setNbLobes(_n);
		
	// Parse the diffuse
	for(int i=0; i<_nY; ++i)
	{
		file >> _kd[i];
	}

	// Parse the lobe
	int n=0;
	std::string line ;
	while(n < _n)
	{
		std::getline(file, line) ;

		if(line.size() > 1)
		{
			std::string sub = line.substr(0,2);

			if(sub == "#C")
			{
				for(int i=0; i<_nY; ++i)
				{
					file >> _C[(n*_nY + i)*2 + 0] >> _C[(n*_nY + i)*2 + 1];
				}
			}
			else if(sub == "#N")
			{
				for(int i=0; i<_nY; ++i)
				{
					file >> _N[n*_nY+i];
				}

				++n;
			}
		}
	}
	
	std::cout << "<<INFO>> Kd = " << _kd << std::endl;
	std::cout << "<<INFO>> Cd = " << _C << std::endl;
	std::cout << "<<INFO>> N = " << _N << std::endl;
}

void isotropic_lafortune_function::save(const std::string& filename) const
{
	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	file << "#DIM " << _nX << " " << _nY << std::endl ;
	file << "#NB_LOBES " << _n << std::endl ;
		
	for(int i=0; i<_nY; ++i)
	{
		file << _kd[i] << std::endl;
	}
	file << std::endl;

	for(int n=0; n<_n; ++n)
	{
		file << "#Lobe number " << n << std::endl;
		file << "#Cxy Cz" << std::endl;
		for(int i=0; i<_nY; ++i)
		{
			file << _C[(n*_nY + i)*2 + 0] << " " << _C[(n*_nY + i)*2 + 1] << std::endl;
		}
		file << std::endl;

		file << "#N" << std::endl;
		for(int i=0; i<_nY; ++i)
		{
			file << _N[n*_nY + i] << std::endl;
		}
		file << std::endl;
	}

}

//! \brief Output the function using a BRDF Explorer formating.
//! \todo Finish
void isotropic_lafortune_function::save_brdfexplorer(const std::string& filename,
                                       const arguments& args) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);
    file << "analytic" << std::endl;

    file << std::endl;
//    file << "::begin parameters" << std::endl;
//    file << "::end parameters" << std::endl;
//    file << std::endl;
    file << std::endl;
    file << "::begin shader" << std::endl;
	 
    type_definition(file, _nY) << " n;" << std::endl;
    type_definition(file, _nY) << " Cx;" << std::endl;
    type_definition(file, _nY) << " Cy;" << std::endl;
    type_definition(file, _nY) << " Cz;" << std::endl;
    type_definition(file, _nY) << " ";
    type_affectation(file, std::string("kd"), _kd, _nY);

    file << std::endl;
    file << std::endl;
    file << "vec3 BRDF( vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y )" << std::endl;
    file << "{" << std::endl;
    if(_nY == 1)
    {
        file << "    ";
        type_definition(file, _nY) << " D = kd;" << std::endl << std::endl;

        for(int n=0; n<_n; ++n)
        {
            file << "    // Lobe number " << n+1 << std::endl;
            file << "    n  = " << _N[n] << "; " << std::endl;
            file << "    Cx = " << _C[2*n + 0] << "; " << std::endl;
            file << "    Cy = " << _C[2*n + 0] << "; " << std::endl;
            file << "    Cz = " << _C[2*n + 1] << "; " << std::endl;
            file << "    D += pow(max(Cx * L.x * V.x + Cy * L.y * V.y + Cz * L.z * V.z, ";
            type_definition(file, _nY) << "(0.0)), n);" << std::endl;
            file << std::endl;
        }
    }
    else
    {
        file << "    ";
        type_definition(file, _nY) << " D = kd;" << std::endl << std::endl;

        for(int n=0; n<_n; ++n)
        {
            file << "    // Lobe number " << n+1 << std::endl;
            file << "    "; type_affectation(file, std::string("n"), _N, _nY, n);
            file << "    "; type_affectation(file, std::string("Cx"), _C, _nY, n, 0, 2);
            file << "    "; type_affectation(file, std::string("Cy"), _C, _nY, n, 0, 2);
            file << "    "; type_affectation(file, std::string("Cz"), _C, _nY, n, 1, 2);
            file << "    D += pow(max(Cx * L.x * V.x + Cy * L.y * V.y + Cz * L.z * V.z, ";
            type_definition(file, _nY) << "(0.0)), n);" << std::endl;
            file << std::endl;
        }
    }
    file << "    return vec3(D);" << std::endl;
//    file << "    if (normalized)" << std::endl;
//    file << "        D *= (2+n) / (2*PI);" << std::endl;
    file << "}" << std::endl;
    file << std::endl;
    file << "::end shader" << std::endl;
    file.close();
}
