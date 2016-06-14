/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/common.h>
#include <core/args.h>

//#define USE_DELAUNAY
#ifdef USE_DELAUNAY
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Delaunay_d.h>


typedef CGAL::Homogeneous_d<double> Kernel;
typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
typedef Delaunay_d::Point_d Point;
typedef Delaunay_d::Simplex_handle S_handle;
typedef Delaunay_d::Point_const_iterator P_iter;

Delaunay_d* D;
int dD;
#else
#include <flann/flann.hpp>
#endif

using namespace alta;

/*! \ingroup datas
 *  \ingroup plugins
 *  \class data_interpolant_rbf
 *  \brief This plugin provide an interpolation of \ref vertical_segment
 *  data using *Radial Basis Function* interpolation.
 *
 *  \details
 *  This interpolation plugin smooth sparse data measurments using radial
 *  kernels. The interpolation is performed as:
 *  <center>
 *     \f$ \tilde{y}(\mathbf{x}) = {\sum_{\mathbf{x}_i \in \mathcal{B}(
 *     \mathbf{x})} k(\mathbf{x} - \mathbf{x}_i) y_i \over \sum_{\mathbf{x}_i
 *     \in \mathcal{B}(\mathbf{x})} k(\mathbf{x} - \mathbf{x}_i)}\f$
 *  </center>
 *
 *  The kernel is set as \f$ k(\mathbf{d}) = {1 \over \epsilon + ||d||^2} \f$.
 *
 *  ### Requirements
 *  This plugin requires the FLANN library to compile. On linux plateforms
 *  it can be obtained by package `libflann-dev` and on OSX using the port
 *  `flann`.
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 */
class rbf_interpolant : public data
{
	private: // data
		// The data object used to load sparse points sets
		ptr<data> _data;

		// Interpolation
#ifndef USE_DELAUNAY
		flann::Index< flann::L2<double> >* _kdtree;
#endif
		int _knn;

	public:

		rbf_interpolant() : _data(new vertical_segment())
		{
			_knn = 3;
		}

		virtual ~rbf_interpolant()
		{
		#ifndef USE_DELAUNAY
			if(_kdtree != NULL)
				delete _kdtree;
		#else
			if(D != NULL)
				delete D;
		#endif
		}

		// Load data from a file
		virtual void load(std::istream& input, const arguments& args)
		{
			// Load the data
      _data->load(input, args);

			// Copy the informations
      parameters p(_data->parametrization().dimX(),
                   _data->parametrization().dimY(),
                   _data->parametrization().input_parametrization(),
                   _data->parametrization().output_parametrization());
      setParametrization(p);
			setMin(_data->min());
			setMax(_data->max());

		#ifdef USE_DELAUNAY
			dD = parametrization().dimX()+parametrization().dimY();
			D  = new Delaunay_d(dD);
			for(int i=0; i<_data->size(); ++i)
			{
				vec x = _data->get(i);

				Point pt(dD, &x[0], &x[dD]);
				D->insert(pt);
			}

			std::cout << "<<DEBUG>> number of points in the Delaunay triangulation: " << D->all_points().size() << std::endl;
			std::cout << "<<DEBUG>> number of points in input: " << _data->size() << std::endl;
		#else
			// Update the KDtreee by inserting all points
			double* _d = new double[parametrization().dimX()*_data->size()];
			flann::Matrix<double> pts(_d, _data->size(), parametrization().dimX());
			for(int i=0; i<_data->size(); ++i)
			{
				vec x = _data->get(i);
				memcpy(pts[i], &x[0], parametrization().dimX()*sizeof(double));
			}
			_kdtree = new flann::Index< flann::L2<double> >(pts, flann::KDTreeIndexParams(4));
			_kdtree->buildIndex();
		#endif
		}

		virtual void save(const std::string& filename) const
		{
		}

		// Acces to data
		virtual vec get(int id) const
		{
			vec res(parametrization().dimX() + parametrization().dimY()) ;
			return res ;
		}
		virtual vec operator[](int i) const
		{
			return get(i) ;
		}

		//! \todo Test this function
		virtual void set(const vec& x)
		{
			NOT_IMPLEMENTED();
		}
		virtual void set(int i, const vec& x)
		{
			NOT_IMPLEMENTED();
		}

		virtual vec value(const vec& x) const
		{
			vec res = vec::Zero(parametrization().dimY());

		#ifndef USE_DELAUNAY
			// Query point
			vec xc(x);
			flann::Matrix<double> pts(&xc[0], 1, parametrization().dimX());
			std::vector< std::vector<int> >    indices;
			std::vector< std::vector<double> > dists;

			_kdtree->knnSearch(pts, indices, dists, _knn, flann::SearchParams());

			// Interpolate the value using the indices
			double cum_dist = 0.0;
			for(int i=0; i<indices[0].size(); ++i)
			{
				const int indice = indices[0][i];
				vec y = _data->get(indice);

		      const double kernel = 1.0/(1.0E-10 + dists[0][i]);

        res      += kernel * y.tail(parametrization().dimY());
				cum_dist += kernel;
			}
			if(cum_dist > 0.0)
			{
				res /= cum_dist;
			}
		#else

			Point pt_x(dD);
			for(int j=0; j<parametrization().dimX(); ++j) { pt_x[j] = x[j]; }
			S_handle simplex = D->locate(pt_x);

			if(simplex == Delaunay_d::Simplex_handle())
			{
				return res;
			}

			// Interpolate the value of the vertices of the Simplex to estimate
			// the value of x
			double cum_dist = 0.0;
			for(int j=0; j<=dD; ++j)
			{
				Point y = D->point_of_simplex(simplex, j);

				// Compute the distance between y and x
				double dist = 0.0;
				for(int i=0; i<parametrization().dimX(); ++i) { dist += pow(y.homogeneous(i)-x[i], 2); }
				dist = sqrt(dist);

				// Interpolate the data
				cum_dist += dist;
				for(int i=0; i<parametrization().dimY(); ++i)
				{
					res[i] += dist * y.homogeneous(parametrization().dimX() + i);
				}

			}

			if(cum_dist > 0.0)
			{
				for(int j=0; j<parametrization().dimY(); ++j)
				{
					res[j] /= cum_dist;
				}
			}
		#endif

		   return res;
		}

		// Get data size, e.g. the number of samples to fit
		virtual int size() const
		{
			assert(_data);
			return _data->size();
		}
};

ALTA_DLL_EXPORT data* provide_data(const arguments&)
{
    return new rbf_interpolant();
}
