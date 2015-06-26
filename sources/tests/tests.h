/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
   Copyright (C) 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/* Utilities for unit tests.  */

#pragma once

#include <cmath>

template<typename T>
static bool close_to(const T a, const T b, const T epsilon = 1E-7)
{
    return std::abs(a - b) < epsilon;
}

template<typename T>
static bool in_range(const T x, const T a, const T b)
{
	return x >= a && x <= b;
}

template<typename T>
static T degrees_to_radians(const T degrees)
{
    static const T d2r = M_PI / 180.;
    return degrees * d2r;
}
