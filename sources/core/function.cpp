#include "function.h"		

#include "common.h"

//! \brief L2 norm to data.
double function::L2_distance(const data* d) const
{
	double l2 = 0.0;
	int nb_points = d->size();
	for(int i=0; i<nb_points; ++i)
	{
		// Get data point
		vec x = d->get(i);
		vec y1(d->dimY());
		for(int j=0; j<d->dimY(); ++j)
			y1 = x[d->dimX() + j];

		// Evaluate data
		vec y2 = this->value(x);

		double dist = norm(y1-y2);
		l2 += dist*dist;
	}

	double factor = 1.0/(double)nb_points;
	vec _min = d->min();
	vec _max = d->max();
	for(int i=0; i<d->dimX(); ++i)
	{
		factor *= _max[i]-_min[i];
	}

	return sqrt(l2)*factor;
}

//! \brief Linf norm to data.
double function::Linf_distance(const data* d) const
{
	double linf = 0.0;
	int nb_points = d->size();
	for(int i=0; i<nb_points; ++i)
	{
				// Get data point
		vec x = d->get(i);
		vec y1(d->dimY());
		for(int j=0; j<d->dimY(); ++j)
			y1 = x[d->dimX() + j];

		// Evaluate data
		vec y2 = this->value(x);

		double dist = norm(y1-y2);
		linf = std::max(dist, linf);
	}

	return linf;
}
