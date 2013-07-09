#include "function.h"		

#include "common.h"

//! \brief L2 norm to data.
double function::L2_distance(const data* d) const
{
	double l2_dist = 0.0;
    for(int i=0; i<d->size(); ++i)
    {
        vec dat = d->get(i);
        vec y(d->dimY());
        for(int j=0; j<d->dimY(); ++j)
            y[j] = dat[d->dimX()+j];

        //linf_dist = std::max<double>(linf_dist, std::abs<double>(norm(y-rj->value(dat))));
        l2_dist  += std::pow(norm(y-value(dat)), 2);
    }
    l2_dist = std::sqrt(l2_dist / d->size());
	return l2_dist;
}

//! \brief Linf norm to data.
double function::Linf_distance(const data* d) const
{

	double linf_dist = 0.0;
    for(int i=0; i<d->size(); ++i)
    {
        vec dat = d->get(i);
        vec y(d->dimY());
        for(int j=0; j<d->dimY(); ++j)
            y[j] = dat[d->dimX()+j];

        linf_dist = std::max<double>(linf_dist, std::abs(norm(y-value(dat))));
    }
	
	return linf_dist;
}
