/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include "common.h"
#include "params.h"

#include <vector>

namespace alta {

#ifdef OLD
class clustering : public data
{
    public: // methods

        //! \brief constructor loading a full dimension data and clustering
        //! it into a low dimension one.
        clustering::clustering(const ptr<data> d, const arguments& args);

        //! \brief the clustering class can save a clustering to a file.
        virtual void save(const std::string& filename);

        //! \brief the clustering class cannot load from a file. It requires
        //! a complete object to be created.
        virtual void load(const std::string& filename)
        {
            throw;
        }
        //! \brief the clustering class cannot load from a file. It requires
        //! a complete object to be created.
        virtual void load(const std::string& filename, const arguments& args)
        {
            throw;
        }

        //! \brief aces to data in linear order
        virtual vec get(int i) const ;
        //! \brief aces to data in linear order
        virtual vec operator[](int i) const ;

        //! \brief return the size of the data after clusterization
        int size() const;

        //! \brief get min input space values
        virtual vec min() const ;
        //! \brief get max input space values
        virtual vec max() const ;

    protected:
        std::vector<vec> _data;
        vec _min, _max;
};
#else
    template<class T> void clustering(const T* in_data , int nY, params::input in_param, params::input out_param, std::vector<vec>& out_data);
#endif

#include "clustering.cpp"
}
