/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

/*! \class ptr
 *
 * Provide an alias for C++11's std::shared_ptr.  Once upon a time this file
 * contained a custom implementation of shared pointers.
 */

#include <memory>

namespace alta {

template<class T> using ptr = std::shared_ptr<T>;

template<class T, class U>
inline ptr<U> dynamic_pointer_cast(const ptr<T>& ptr_t)
{
    return std::dynamic_pointer_cast<U>(ptr_t);
}
}
