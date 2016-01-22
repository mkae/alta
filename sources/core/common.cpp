/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "common.h"

#ifdef _WIN32
#include <windows.h>
#include <winbase.h>
#else
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>
#include <unistd.h>
#endif

#include <iostream>
#include <cassert>

using namespace alta;

vec product(const vec& a, const vec& b)
{
    if (a.size() == 1)
        return b * a[0];
    else if (b.size() == 1)
        return a * b[0];
    else
        return a.cwiseProduct(b);
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
    _elapsed += _stop - _start;
    _start = _stop;
}

void timer::reset()
{
    _elapsed = 0;
    _start   = current_time();
    _stop    = _start;
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
#if defined(_WIN32)
	SYSTEMTIME res;
	GetSystemTime(&res);
	return (unsigned int)(res.wSecond + res.wMinute*60 + res.wHour*360);
#elif defined(__APPLE__)
    struct timeval _time;
    gettimeofday(&_time, NULL);
    return (unsigned int)_time.tv_sec;
#else	
	struct timespec _time;
	clock_gettime(CLOCK_MONOTONIC, &_time);
	return (unsigned int)_time.tv_sec;
#endif
}

