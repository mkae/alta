#include "common.h"

double norm(const vec& a)
{
	double norm = 0.0 ;
	for(int i=0; i<a.size(); ++i)
	{
		norm += a[i]*a[i];
	}
	return sqrt(norm);
}

vec normalize(const vec& a)
{
	vec b(a.size());
	double norm = 0.0 ;
	for(int i=0; i<a.size(); ++i)
	{
		norm += a[i]*a[i];
	}
	norm = sqrt(norm);

	for(int i=0; i<a.size(); ++i)
	{
		b[i] = a[i]/norm;
	}
	return b;
}

double dot(const vec& a, const vec& b)
{
#ifdef DEBUG
	assert(a.size() == b.size());
#endif
	double res = 0.0;
	for(int i=0; i<a.size(); ++i)
	{
		res += a[i]*b[i];
	}

	return res;
}

vec product(const vec& a, const vec& b)
{
#ifdef DEBUG
	assert(a.size() == b.size());
#endif
	vec res(a.size());

	for(int i=0; i<a.size(); ++i)
	{
		res[i] = a[i]*b[i];
	}

	return res;
}

std::ostream& operator<<(std::ostream& out, const vec& v)
{
    out << "[";
    for(int i=0; i<v.size(); ++i)
    {
        out << v(i);
        if(i != v.size()-1) { out << ", "; }
    }
    out << "]";
    return out;
}

