/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
   Copyright (C) 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/* Utilities for unit tests.  */

#pragma once

// Make sure 'assert' does its job in tests.
#undef NDEBUG

#include <iostream>
#include <cstdlib>
#include <cmath>

namespace alta { namespace tests
{

#define STRINGIFY_(x) #x
#define STRINGIFY(x)  STRINGIFY_(x)

#define TEST_ASSERT(exp)                                          \
    do																														\
    {																															\
        std::cerr << "evaluating '" << STRINGIFY(exp) << "'... ";	\
        bool result = (exp);																			\
        std::cerr << (result ? "PASS" : "FAIL") << std::endl;			\
        if (!result)																							\
            abort();																							\
    }																															\
    while(0)

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

    template<typename T>
    class angle_range
    {
    public:
        class angle_iterator
        {
        public:
            angle_iterator(const T start, const T step):
            _value(start), _step(step) { }

            bool operator!=(const angle_iterator& other) const
            {
                return _value != other._value;
            }

            T operator*() const
            {
                return degrees_to_radians(_value);
            }

            const angle_iterator& operator++()
            {
                _value += _step;
                return *this;
            }

        private:
            T _value;
            const T _step;
        };

        angle_range(T start, T end, T step)
            : _start(start), _end(end), _step(step)
        { }

        angle_iterator begin()
        {
            return angle_iterator(_start, _step);
        }

        angle_iterator end()
        {
            return angle_iterator(_end + _step, _step);
        }

    private:
        const T _start, _end, _step;
    };
}};
