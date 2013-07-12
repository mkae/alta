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
    return new lafortune_function();
}

// Overload the function operator
vec lafortune_function::operator()(const vec& x) const 
{
	return value(x);
}
vec lafortune_function::value(const vec& x) const 
{
	vec res(dimY());

#ifdef ADAPT_TO_PARAM
	vec y(6);
	params::convert(&x[0], _in_param, params::CARTESIAN, &y[0]);
#endif

	double dx, dy, dz;
#ifdef ADAPT_TO_PARAM
	dx = y[0]*y[3];
	dy = y[1]*y[4];
	dz = y[2]*y[5];
#else
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];
#endif

	// For each color channel
	for(int i=0; i<dimY(); ++i)
	{
		// Start with the diffuse term
		res[i] = _kd[i];

		// For each lobe
		for(int n=0; n<_n; ++n)
		{
			double Cx, Cy, Cz, N;
			getCurrentLobe(n, i, Cx, Cy, Cz, N);

			const double d = Cx*dx + Cy*dy + Cz*dz;
			if(d > 0.0)
				res[i] += pow(d, N);
		}

	}

    return res;
}
        
vec lafortune_function::value(const vec& x, const vec& p) const
{
	// Test input parameters for correct size
	assert(p.size() == nbParameters());

	vec res(dimY());

#ifdef ADAPT_TO_PARAM
	vec y(6);
	params::convert(&x[0], _in_param, params::CARTESIAN, &y[0]);
#endif

	double dx, dy, dz;
#ifdef ADAPT_TO_PARAM
	dx = y[0]*y[3];
	dy = y[1]*y[4];
	dz = y[2]*y[5];
#else
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];
#endif

	// For each lobe and for each color channel
	for(int i=0; i<dimY(); ++i)
	{
		// Start with the diffuse
		res[i] = _kd[i];

		// Add the n lobes
		for(int n=0; n<_n; ++n)
		{
			double N, Cx, Cy, Cz;
			getCurrentLobe(n, i, Cx, Cy, Cz, N);

			const double d = Cx*dx + Cy*dy + Cz*dz;
			if(d > 0.0)
				res[i] += pow(d, N);
		}
	}

	return res;
}

// Set the number of lobes of the Lafortune BRDF
void lafortune_function::setNbLobes(int N)
{
    _n = N;

    // Update the length of the vectors
    _C.resize(_n*_nY*3) ;
    _N.resize(_n*_nY) ;
}

// Reset the output dimension
void lafortune_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _C.resize(_n*_nY*3) ;
    _N.resize(_n*_nY) ;
    _kd.resize(_nY);

    for(int i=0; i<nY; ++i)
        _kd[i] = 0.0;
}

//! Number of parameters to this non-linear function
int lafortune_function::nbParameters() const 
{
#ifdef FIT_DIFFUSE
    return (4*_n+1)*dimY();
#else
    return (4*_n)*dimY();
#endif
}

//! Get the vector of parameters for the function
vec lafortune_function::parameters() const 
{
#ifdef FIT_DIFFUSE
    vec res((4*_n+1)*dimY());
#else
    vec res((4*_n)*dimY());
#endif
    for(int n=0; n<_n; ++n)
	    for(int i=0; i<dimY(); ++i)
		 {
			  res[(n*dimY() + i)*4 + 0] = _C[(n*dimY() + i)*3 + 0];
			  res[(n*dimY() + i)*4 + 1] = _C[(n*dimY() + i)*3 + 1];
			  res[(n*dimY() + i)*4 + 2] = _C[(n*dimY() + i)*3 + 2];
			  res[(n*dimY() + i)*4 + 3] = _N[n*dimY()  + i];
		 }

#ifdef FIT_DIFFUSE
	 for(int i=0; i<dimY(); ++i)
	 {
		 res[4*_n*dimY() + i] = _kd[i];
	 }
#endif
    return res;
}

//! Update the vector of parameters for the function
void lafortune_function::setParameters(const vec& p) 
{
	// Safety check the number of parameters
	assert(p.size() == nbParameters());

	for(int n=0; n<_n; ++n)
		for(int i=0; i<dimY(); ++i)
		{
			_C[(n*dimY() + i)*3 + 0] = p[(n*dimY() + i)*4 + 0];
			_C[(n*dimY() + i)*3 + 1] = p[(n*dimY() + i)*4 + 1];
			_C[(n*dimY() + i)*3 + 2] = p[(n*dimY() + i)*4 + 2];
			_N[n*dimY()  + i]        = p[(n*dimY() + i)*4 + 3];
		}
#ifdef FIT_DIFFUSE
	for(int i=0; i<dimY(); ++i)
	{
		_kd[i] = p[4*_n*dimY() + i];
	}
#endif
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec lafortune_function::parametersJacobian(const vec& x) const 
{

#ifdef ADAPT_TO_PARAM
	vec y(6);
	params::convert(&x[0], _in_param, params::CARTESIAN, &y[0]);
#endif

	double dx, dy, dz;
#ifdef ADAPT_TO_PARAM
	dx = y[0]*y[3];
	dy = y[1]*y[4];
	dz = y[2]*y[5];
#else
	dx = x[0]*x[3];
	dy = x[1]*x[4];
	dz = x[2]*x[5];
#endif

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int n=0; n<_n; ++n)
			 for(int j=0; j<dimY(); ++j)
			 {
				 // index of the current monochromatic lobe
				 int index = i*nbParameters() + 4*(n*dimY() + j);
				 
				 double Cx, Cy, Cz, N;
				 getCurrentLobe(n, j, Cx, Cy, Cz, N);
				 
				 double d  = Cx*dx + Cy*dy + Cz*dz;

				 if(i == j && d > 0.0)
				 {
					 // df / dCx
					 jac[index+0] = dx * N * std::pow(d, N-1.0);

					 // df / dCy
					 jac[index+1] = dy * N * std::pow(d, N-1.0);
					 
					 // df / dCz
					 jac[index+2] = dz * N * std::pow(d, N-1.0);

					 // df / dN
					 if(d <= 0.0)
						 jac[index+3] = 0.0;
					 else
						 jac[index+3] = std::log(d) * std::pow(d, N);
				 }
				 else
				 {
					 jac[index+0] = 0.0;
					 jac[index+1] = 0.0;
					 jac[index+2] = 0.0;
					 jac[index+3] = 0.0;
				 }
			 }

#ifdef FIT_DIFFUSE
		 for(int j=0; j<dimY(); ++j)
		 {
			 // index of the current monochromatic lobe
			 int index = i*nbParameters() + 4*_n*dimY() + j;

			 jac[index] = 1.0;
		 }
#endif
	 }

    return jac;
}
		
void lafortune_function::bootstrap(const data* d, const arguments& args)
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

    // The remaining data will be equal to one
    for(int n=0; n<_n; ++n)
        for(int i=0; i<dimY(); ++i)
        {
            double theta = 0.5 * M_PI * n / (double)_n;

            _C[(n*dimY() + i)*3 + 0] = -sin(theta);
            _C[(n*dimY() + i)*3 + 1] = -sin(theta);
            _C[(n*dimY() + i)*3 + 2] = cos(theta);
            _N[n*dimY()  + i]        = (double)_n;
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
void lafortune_function::load(const std::string& filename)
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

	setNbLobes(_n);

	for(int n=0; n<_n; ++n)
	{
		std::cout << (char)file.peek() << std::endl;

		for(int i=0; i<_nY; ++i)
		{
			file >> _C[(n*_nY + i)*3 + 0] >> _C[(n*_nY + i)*3 + 1] >> _C[(n*_nY + i)*3 + 2];
		}
	}
}

void lafortune_function::save(const std::string& filename) const
{
	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	file << "#DIM " << _nX << " " << _nY << std::endl ;
	file << "#NB_LOBES " << _n << std::endl ;

	for(int n=0; n<_n; ++n)
	{
		file << "#Lobe number " << n << std::endl;
		file << "#Cx Cy Cz" << std::endl;
		for(int i=0; i<_nY; ++i)
		{
			file << _C[(n*_nY + i)*3 + 0] << " " << _C[(n*_nY + i)*3 + 2] << " " << _C[(n*_nY + i)*3 + 1] << std::endl;
		}
		file << std::endl;

		file << "#N" << std::endl;
		for(int i=0; i<_nY; ++i)
		{
			file << _N[i] << std::endl;
		}
		file << std::endl;
	}

}

//! \brief Output the function using a BRDF Explorer formating.
//! \todo Finish
void lafortune_function::save_brdfexplorer(const std::string& filename,
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
            file << "    Cx = " << _C[0*_n + n] << "; " << std::endl;
            file << "    Cy = " << _C[1*_n + n] << "; " << std::endl;
            file << "    Cz = " << _C[2*_n + n] << "; " << std::endl;
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
            file << "    "; type_affectation(file, std::string("Cx"), _C, _nY, n, 0, 3);
            file << "    "; type_affectation(file, std::string("Cy"), _C, _nY, n, 1, 3);
            file << "    "; type_affectation(file, std::string("Cz"), _C, _nY, n, 2, 3);
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


Q_EXPORT_PLUGIN2(lafortune_function, lafortune_function)
