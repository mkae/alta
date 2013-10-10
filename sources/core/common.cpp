#include "common.h"

#ifdef _WIN32
#include <windows.h>
#include <winbase.h>
#else
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>
#include <unistd.h>
#endif

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
    if(a.size() == 1 && b.size() > 1)
    {
        vec res(b.size());
        for(int i=0; i<b.size(); ++i)
        {
            res[i] = a[0]*b[i];
        }
        return res;
    }
    else if(b.size() == 1 && a.size() > 1)
    {
        vec res(a.size());
        for(int i=0; i<a.size(); ++i)
        {
            res[i] = a[i]*b[0];
        }
        return res;
    }
    else
    {
#ifndef DEBUG
        assert(a.size() == b.size());
#endif
        vec res(b.size());
        for(int i=0; i<b.size(); ++i)
        {
            res[i] = a[i]*b[i];
        }
        return res;
    }
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








/* ---- Timer implementation ---- */

timer::timer()
{
    reset();
}

void timer::start()
{
    _elapsed += _stop-_start;
    _start = current_time();
}

void timer::stop()
{
    _stop = current_time();
}

void timer::reset()
{
    _elapsed = 0;
    _start   = 0;
    _stop    = 0;
}

int timer::elapsed() const
{
    return _elapsed;
}

void timer::print(std::ostream& out) const
{
    unsigned int tsec = elapsed() ;
    unsigned int sec  = tsec % 60 ;
    unsigned int min  = (tsec / 60) % 60 ;
    unsigned int hour = (tsec / 360) ;
    out << hour << "h " << min << "m " << sec << "s" ;
}

std::ostream& operator<<(std::ostream& out, const timer& t)
{
    t.print(out);
    return out;
}

unsigned int timer::current_time() const
{
#ifdef WIN32
	SYSTEMTIME res;
	GetSystemTime(&res);
	return (unsigned int)(res.wSecond + res.wMinute*60 + res.wHour*360);
#else
    struct timespec res;
    clock_gettime(CLOCK_MONOTONIC, &res);
    return static_cast<unsigned int>(res.tv_sec);
#endif
}
