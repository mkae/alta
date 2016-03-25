/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <iostream>
#include <cmath>
#include "data.h"
#include "data_storage.h"

using namespace alta;

void data::save(const std::string& filename) const
{
    std::ofstream file;

    file.exceptions(std::ios_base::failbit);
    file.open(filename.c_str(), std::ios_base::trunc);
    file.exceptions(std::ios_base::goodbit);

        alta::save_data_as_text(file, *this);

    file.close();
}

bool data::equals(const data& data, double epsilon)
{
    if (size() != data.size()
        || _parameters.dimX() != data.parametrization().dimX()
        || _parameters.dimY() != data.parametrization().dimY()
        || (_parameters.input_parametrization()
            != data.parametrization().input_parametrization())
        || (_parameters.output_parametrization()
            != data.parametrization().output_parametrization()))
        return false;

    for(int i = 0; i < data.size(); i++)
    {
        vec other = data.get(i);
        vec self = get(i);

        if (self.size() != other.size())          // should not happen
            return false;

        for (int j = 0; j < self.size(); j++)
        {
            double diff = std::abs(self[j] - other[j]);
            if (diff > epsilon)
                return false;
        }
    }

    return true;
}

void data::setMin(const vec& min)
{
#ifdef DEBUG
    assert(min.size() == _nX) ;
#endif
    _min = min ;
}
void data::setMax(const vec& max)
{
#ifdef DEBUG
    assert(max.size() == _nX) ;
#endif
    _max = max ;
}
