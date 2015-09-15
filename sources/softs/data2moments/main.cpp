/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/*! \package data2moment
 *  \ingroup commands
 *  \brief
 *  This command computes the moments of a \ref data object.
 *  \details
 *  <h3>Parameters</h3>
 *  <ul>
 *		<li><b>\-\-input <i>filename</i></b> data file to be loaded uwing the data
 *		plugin specified by the <b>\-\-data <i>filename</i></b> option.
 *		</li>
 *		<li><b>\-\-data <i>filename</i></b> specify the data plugin used to load
 *		the data file. \note It is required to provide a data plugin that performs
 *		interpolation of the data.
 *		</li>
 *  </ul>
 */
#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>
#include <core/vertical_segment.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;

	if(args.is_defined("help")) {
		std::cout << "Usage: data2moments [options] --input data.file --output gnuplot.file" << std::endl ;
		std::cout << "Compute the statistical moments, up to order 4, on the data file" << std::endl ;
		std::cout << std::endl;
		std::cout << "Mandatory arguments:" << std::endl;
		std::cout << "  --input   [filename]" << std::endl;
		std::cout << "  --output  [filename]" << std::endl;
		std::cout << "  --data    [filename]       You must provide a data loader that provides" << std::endl;
		std::cout << "                             interpolation mechanism." << std::endl;
		std::cout << std::endl;
		std::cout << "Optional arguments:" << std::endl;
		std::cout << "  --samples [int, int, ...]  Number of samples per dimension used to compute" << std::endl;
		std::cout << "                             the moments. This vector must have the same" << std::endl;
		std::cout << "                             size as the number of dimensions of the input" << std::endl;
		std::cout << "                             space." << std::endl;
		std::cout << "  --dim     [int, int]       Indices of the two dimensions used to evaluate" << std::endl;
		std::cout << "                             the moments." << std::endl;
		return 0 ;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("data")) {
		std::cerr << "<<ERROR>> the data loader is not defined" << std::endl ;
		return 1 ;
	}

	// Import data
	ptr<data> d = NULL ;
	d = plugins_manager::get_data(args["data"], args) ;

	if(dynamic_pointer_cast<vertical_segment>(d)) {
		std::cerr << "<<ERROR>> this data object is not interpolant." << std::endl;
		return 1;
	}
	try
	{
		d->load(args["input"]);
	}
	CATCH_FILE_IO_ERROR(args["input"]);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	if(d)
	{
		// Data parametrization
		/*params::input data_param = d->parametrization();*/
		const int nX = d->dimX();
		const int nY = d->dimY();

		// Raw moments
		vec m_0(nY);
		vec m_x(nY);
		vec m_y(nY);
		vec m_xx(nY);
		vec m_xy(nY);
		vec m_yy(nY);
		vec m_xxx(nY);
		vec m_xxy(nY);
		vec m_xyy(nY);
		vec m_yyy(nY);
		vec m_xxxx(nY);
		vec m_xxxy(nY);
		vec m_xxyy(nY);
		vec m_xyyy(nY);
		vec m_yyyy(nY);

		// Cumulants
		vec k_x(nY);
		vec k_y(nY);
		vec k_xx(nY);
		vec k_xy(nY);
		vec k_yy(nY);
		vec k_xxx(nY);
		vec k_xxy(nY);
		vec k_xyy(nY);
		vec k_yyy(nY);
		vec k_xxxx(nY);
		vec k_xxxy(nY);
		vec k_xxyy(nY);
		vec k_xyyy(nY);
		vec k_yyyy(nY);

		// Number of elements per samplesension
		std::vector<int> samples;
		if(args.is_defined("samples"))
		{	
			samples = args.get_vec<int>("samples");
			assert(samples.size() == nX);
		}
		else
		{
			for(int i=0; i<nX; ++i)
			{
				samples.push_back(100);
			}
		}

		std::vector<int> dim;
		if(args.is_defined("dim"))
		{
			dim = args.get_vec<int>("dim");
			assert(dim.size() == 2);
		}
		else
		{
			std::cout << "<<INFO>> not dim [int, int] defined. Will use the first two dimensions to compute the moments" << std::endl;
			dim.push_back(0);
			dim.push_back(1);
		}


		// Get the minimum element of the dataset
		vec _d_min = d->get(0);
		for(int i=1; i<d->size(); ++i)
		{
			const vec y = d->get(i);
			for(int j=0; j<nY; ++j)
			{
				_d_min[nX + j] = std::min(_d_min[nX + j], y[nX + j]);
			}
		}

		// Compute the volume in which we integrate and then compute the
		// dt to apply for each element.
		double dt = 1.0;
		int nb = 1;
		for(int k=0; k<nX; ++k)
		{
			const double length = d->max()[k] - d->min()[k];

			if(length > 0.0)
			{
				dt *= length;
				nb *= samples[k];
			}
			else
			{
				samples[k] = 1;
			}
		}
		std::cout << "<<DEBUG>> volume of the integration domain: " << dt << std::endl;
		dt /= double(nb);
		std::cout << "<<DEBUG>> volume of an integration step: " << dt << std::endl;
		
		// Set all values to zero
		m_0  = vec::Zero(nY);
		m_x  = vec::Zero(nY);
		m_y  = vec::Zero(nY);
		m_xx = vec::Zero(nY);
		m_xy = vec::Zero(nY);
		m_yy = vec::Zero(nY);
		m_xxx = vec::Zero(nY);
		m_xxy = vec::Zero(nY);
		m_xyy = vec::Zero(nY);
		m_yyy = vec::Zero(nY);
		m_xxxx = vec::Zero(nY);
		m_xxxy = vec::Zero(nY);
		m_xxyy = vec::Zero(nY);
		m_xyyy = vec::Zero(nY);
		m_yyyy = vec::Zero(nY);

		const int dx = dim[0];
		const int dy = dim[1];
		int* indices = new int[nX];

		// Integrate the moments over the integration domain
		for(int i=0; i<nb; ++i)
		{
			// For the global index i of the sample, compute the index in the various
			// dimensions and then compute the position within the dimension. All the
			// samples are taken uniformly.
			int global = i;
			for(int k=0; k<nX; ++k)
			{
				indices[k] = global % samples[k];
				global /= samples[k];
			}

			vec x(nX);
			for(int k=0; k<nX; ++k)
			{
				x[k] = d->min()[k] + (d->max()[k] - d->min()[k]) * (double(indices[k]) / samples[k]);
			}

			// Get the value and compute the associated integral. Right now, the data
			// is supposed to be 2-dimensional.
			// \todo Handle multi-dimensional data
			vec y = d->value(x);
			for(int k=0; k<nY; ++k)
			{
				y[k] -= _d_min[nX + k];

				double val = y[k] * dt;

				m_0[k]  += val ;
				m_x[k]  += val * x[dx];
				m_y[k]  += val * x[dy];
				m_xx[k] += val * x[dx] * x[dx];
				m_xy[k] += val * x[dx] * x[dy];
				m_yy[k] += val * x[dy] * x[dy];
				m_xxx[k] += val* x[dx]  * x[dx] * x[dx];
				m_xxy[k] += val* x[dx]  * x[dx] * x[dy];
				m_xyy[k] += val* x[dx]  * x[dy] * x[dy];
				m_yyy[k] += val* x[dy]  * x[dy] * x[dy];
				m_xxxx[k] += val* x[dx]  * x[dx] * x[dx] * x[dx];
				m_xxxy[k] += val* x[dx]  * x[dx] * x[dx] * x[dy];
				m_xxyy[k] += val* x[dx]  * x[dx] * x[dy] * x[dy];
				m_xyyy[k] += val* x[dx]  * x[dy] * x[dy] * x[dy];
				m_yyyy[k] += val* x[dy]  * x[dy] * x[dy] * x[dy];
			}
		}
		
		// Clean memory
		delete[] indices;

		for(int k=0; k<nY; ++k)
		{
			m_x[k]  /= m_0[k];
			m_y[k]  /= m_0[k];
			m_xx[k] /= m_0[k];
			m_xy[k] /= m_0[k];
			m_yy[k] /= m_0[k];
			m_xxx[k] /= m_0[k];
			m_xxy[k] /= m_0[k];
			m_xyy[k] /= m_0[k];
			m_yyy[k] /= m_0[k];
			m_xxxx[k] /= m_0[k];
			m_xxxy[k] /= m_0[k];
			m_xxyy[k] /= m_0[k];
			m_xyyy[k] /= m_0[k];
			m_yyyy[k] /= m_0[k];
		}

		// compute cumulants
		k_x  = m_x;
		k_y  = m_y;
		k_xx = m_xx - product(m_x,m_x);
		k_xy = m_xy - product(m_x,m_y);
		k_yy = m_yy - product(m_y,m_y);
		k_xxx = m_xxx - 3*product(m_xx, m_x) + 2*product(m_x, product(m_x, m_x));
		k_xxy = m_xxy - product(m_xx, m_y) - 2*product(m_xy, m_x) +2*product(m_x, product(m_x, m_y));
		k_xyy = m_xyy - product(m_yy, m_x) - 2*product(m_xy, m_y) +2*product(m_x, product(m_y, m_y));
		k_yyy = m_yyy - 3*product(m_yy, m_y)+2*product(m_y, product(m_y, m_y));
#ifdef NOT
		k_xxxx = m_xxxx - 4*m_xxx*m_x - 3*m_xx*m_xx + 12*m_xx*m_x*m_x - 6*m_x*m_x*m_x*m_x;
		k_xxxy = m_xxxy - 3*m_xxy*m_x - m_xxx*m_y   - 3*m_xx*m_xy     + 6*m_xx*m_x*m_y + 6*m_xy*m_x*m_x - 6*m_x*m_x*m_x*m_y;
		k_xxyy = m_xxyy - 2*m_xxy*m_y - 2*m_xyy*m_x   - 2*m_xy*m_xy     - m_xx*m_yy      + 8*m_xy*m_x*m_y + 2*m_xx*m_y*m_y + 2*m_yy*m_x*m_x - 6*m_x*m_x*m_y*m_y;
		k_xyyy = m_xyyy - 3*m_xyy*m_y - m_yyy*m_x   - 3*m_yy*m_xy     + 6*m_yy*m_x*m_y + 6*m_xy*m_y*m_y - 6*m_x*m_y*m_y*m_y;
		k_yyyy = m_yyyy - 4*m_yyy*m_y - 3*m_yy*m_yy + 12*m_yy*m_y*m_y - 6*m_y*m_y*m_y*m_y;
#endif
		std::cout << "<<RESULT>> moment    0: " << m_0  << std::endl;

		std::cout << "<<RESULT>> moment    x: " << m_x  << std::endl;
		std::cout << "<<RESULT>> moment    y: " << m_y  << std::endl;

		std::cout << "<<RESULT>> moment   xx: " << m_xx << std::endl;
		std::cout << "<<RESULT>> moment   xy: " << m_xy << std::endl;
		std::cout << "<<RESULT>> moment   yy: " << m_yy << std::endl;

		std::cout << "<<RESULT>> moment  xxx: " << m_xxx << std::endl;
		std::cout << "<<RESULT>> moment  xxy: " << m_xxy << std::endl;
		std::cout << "<<RESULT>> moment  xyy: " << m_xyy << std::endl;
		std::cout << "<<RESULT>> moment  yyy: " << m_yyy << std::endl;
#ifdef NOT
		std::cout << "<<RESULT>> moment xxxx: " << m_xxxx << std::endl;
		std::cout << "<<RESULT>> moment xxxy: " << m_xxxy << std::endl;
		std::cout << "<<RESULT>> moment xxyy: " << m_xxyy << std::endl;
		std::cout << "<<RESULT>> moment xyyy: " << m_xyyy << std::endl;
		std::cout << "<<RESULT>> moment yyyy: " << m_yyyy << std::endl;
#endif

		file << "\\mu_0 = " << m_0 << std::endl;
		file << "\\mu_x = "  << m_x << std::endl;
		file << "\\mu_y = "  << m_y << std::endl;
		file << "\\mu_{xx} = "  << m_xx << std::endl;
		file << "\\mu_{xy} = "  << m_xy << std::endl;
		file << "\\mu_{yy} = "  << m_yy << std::endl;
		file << "\\mu_{xxx} = "  << m_xxx << std::endl;
		file << "\\mu_{xxy} = "  << m_xxy << std::endl;
		file << "\\mu_{xyy} = "  << m_xyy << std::endl;
		file << "\\mu_{yyy} = "  << m_yyy << std::endl;
#ifdef NOT
		file << "\\mu_{xxxx} = "  << m_xxxx << std::endl;
		file << "\\mu_{xxxy} = "  << m_xxxy << std::endl;
		file << "\\mu_{xxyy} = "  << m_xxyy << std::endl;
		file << "\\mu_{xyyy} = "  << m_xyyy << std::endl;
		file << "\\mu_{yyyy} = "  << m_yyyy << std::endl;
#endif

		file.close();
		return 0 ;
	}
	else
	{
		file.close();
		std::cerr << "<<ERROR>> unable to load file \"" << args["input"] << "\"" << std::endl ;
		return 1;
	}
}
