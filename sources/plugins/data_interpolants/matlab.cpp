/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2016 Inria
   Copyright (C) 2015 CNRS

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// STL includes
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

// ALTA includes
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/common.h>
#include <core/args.h>

// Matlab includes
#include <engine.h>

using namespace alta;

mxArray *X, *Y, *x, *y;

#define BUFFER_SIZE 10000
char* output = new char[BUFFER_SIZE+1];


/*! \ingroup datas
 *  \ingroup plugins
 *  \class data_interpolant_matlab
 *  \brief This plugin provide an interpolation of \ref vertical_segment
 *  data using Matlab `griddata` functions.
 *
 *  \details
 *  This plugin provide an interpolation of \ref vertical_segment
 *  data using Matlab `griddata` functions.
 *
 *  ### Requirements
 *  This plugin requires the Matlab engine library to compile.
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 *  \author Romain Pacanowski <romain.pacanowski@institutoptique.fr>
 *
 */
class MatlabInterpolant : public data
{
	private: // data

		// The data object used to load sparse points sets
		ptr<data> _data;
		Engine *ep;

	public: // methods
		MatlabInterpolant()
    : _data( ptr<data>( new vertical_segment() ) )
		{
			// Create matlab engine
		#ifdef WIN32
		    if (!(ep = engOpen(NULL)))
		#else
		    if (!(ep = engOpen("matlab -nosplash")))
		#endif
			{
				std::cerr << "<ERROR>> can't start MATLAB engine" << std::endl ;
			}

			output[BUFFER_SIZE] = '\0';
			engOutputBuffer(ep, output, BUFFER_SIZE) ;
		}

    virtual ~MatlabInterpolant()
		{
			delete[] output;

			mxDestroyArray(x);
			mxDestroyArray(y);
			mxDestroyArray(X);
			mxDestroyArray(Y);
			engClose(ep);
		}

		// Load data from a file
		virtual void load(std::istream& input, const arguments& args)
		{
			// Load the data
      _data->load(input, args);

			// Copy the informations
			setDimX(_data->dimX());
			setDimY(_data->dimY());
			setMin(_data->min());
			setMax(_data->max());

			Eigen::MatrixXd eX(_data->size(), dimX()*dimY());
			Eigen::MatrixXd eY(_data->size(), dimY());
			for(int i=0; i<_data->size(); ++i)
			{
				vec x = _data->get(i);

				// For a point extract all its coordinates
				for(int j=0; j<dimX(); ++j)
				{
					eX(i, j) =  x[j];
				}

				// Extract all its values
				for(int k=0; k<dimY(); ++k)
				{
					eY(i, k) = x[dimX() + k];
				}
			}

			// Create matrices
			X = mxCreateDoubleMatrix(_data->size(), dimX(), mxREAL);
			memcpy((void *)mxGetPr(X), (void *) eX.data(),  _data->size()*dimY()*dimX()*sizeof(double));
			engPutVariable(ep, "X", X);

			Y = mxCreateDoubleMatrix(_data->size(), dimY(), mxREAL);
			memcpy((void *)mxGetPr(Y), (void *) eY.data(),  _data->size()*dimY()*sizeof(double));
			engPutVariable(ep, "Y", Y);

			x = mxCreateDoubleMatrix(1, dimX(), mxREAL);
		}

		virtual void save(const std::string& filename) const
		{
			_data->save(filename);
		}

		// Acces to data
		virtual vec get(int id) const
		{
			return _data->get(id);
		}

    inline virtual vec operator[](int i) const
		{
			return get(i) ;
		}

		virtual void set(const vec& x)
		{
			NOT_IMPLEMENTED();
		}
		virtual void set(int i, const vec& x)
		{
			_data->set(i, x);
		}

		virtual vec value(const vec& ax) const
		{
			vec res = vec::Zero(dimY());

			// Copy the input vector to matlab code
			memcpy((void *)mxGetPr(x), (void *)&ax[0],  dimX()*sizeof(double));
			engPutVariable(ep, "x", x);

			std::stringstream cmd;

			// Iterate over the output dimension
			for(int i=0; i<dimY(); ++i)
			{
				// Evaluate the matlab routine
				if(dimX() == 2)
				{
					cmd << "y(" << i+1 << ") = griddata(X(:,1), X(:,2), Y(:," << i+1 << "), x(1), x(2), 'cubic');";
				}
				else if(dimX() == 3)
				{
					cmd << "y(" << i+1 << ") = griddata(X(:,1), X(:,2), X(:,3), Y(:," << i+1 << "), x(1), x(2), x(3), 'cubic');";
				}
				else
				{
					cmd << "y(" << i+1 << ") = griddatan(X(:, 1:" << dimX() <<"), Y(:, " << i+1 << "), x, 'linear');";
				}
				engEvalString(ep, cmd.str().c_str());
		#ifdef DEBUG
				std::cout << output;
		#endif
			}

			// Get results and copy it
			y = engGetVariable(ep, "y") ;
			double* y_val = (double*)mxGetData(y);
			memcpy(&res[0], y_val, dimY()*sizeof(double));

			// Fail safe: if the query point is outside of the convex hull of the
			// data, rerun the algorithm using a nearest option.
			if(isnan(res[0]))
			{
				for(int i=0; i<dimY(); ++i)
				{
					cmd.str(std::string());
					cmd << "y("<< i+1 <<") = griddatan(X(:, 1:" << dimX() <<"), Y(:, " << i+1 << "), x, 'nearest');";
					engEvalString(ep, cmd.str().c_str());
				}
				// Get results and copy it
				y = engGetVariable(ep, "y") ;
				double* y_val = (double*)mxGetData(y);
				memcpy(&res[0], y_val, dimY()*sizeof(double));
			}

		   return res;
		}

		// Get data size, e.g. the number of samples to fit
		virtual int size() const
		{
			return _data->size();
		}
};

ALTA_DLL_EXPORT data* provide_data(const arguments&)
{
    return new MatlabInterpolant();
}
