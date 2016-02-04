/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
   Copyright (C) 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// ALTA includes
#include <core/params.h>
#include <tests.h>

using namespace alta::tests;

// STL includes
#include <cmath>
#include <iostream>
#include <functional>


// Return true if every element of V matches PRED, false otherwise.
static bool every(const std::function<bool (double)> &pred, const vec &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        if (!pred(v[i]))
            return false;
    }
    return true;
}

/* Test different configurations for the Half / Cartesian parametrization
 * Returns: 0 is every test passes
 *          n > 0 when n tests did no pass
 */
int main(int argc, char** argv) {

  // Number of failed tests
  int n = 0;

  const int nphi   = 360;
  const int ntheta = 90;
  const int step   = 10;

  // Sample the full hemisphere for both ωi and ωo.  Check if the converted
  // vector is within the range of the parametrization.  θh and θd should be
  // within [0..π/2], and Φd should be within [-π..π].
  for (auto&& theta_i : angle_range<double>(0, ntheta, step)) {
  for (auto&& phi_i : angle_range<double>(0, nphi, step)) {

      for (auto&& theta_o : angle_range<double>(0, ntheta, step)) {
      for (auto&& phi_o : angle_range<double>(0, nphi, step)) {

          vec cart(6);
          cart[0] = cos(phi_i)*sin(theta_i);
          cart[1] = sin(phi_i)*sin(theta_i);
          cart[2] = cos(theta_i);
          cart[3] = cos(phi_o)*sin(theta_o);
          cart[4] = sin(phi_o)*sin(theta_o);
          cart[5] = cos(theta_o);

          // Are all the elements of CART positive?  IOW, check whether the
          // input coordinates are in the hemisphere.
          if (!every([](double x) { return x >= 0; }, cart))
              // No they aren't: skip this value.
              continue;

          vec x(4);
          params::convert(&cart[0], params::CARTESIAN, params::RUSIN_TH_PH_TD_PD, &x[0]);

          if(!in_range(x[0], 0.0, M_PI_2) || !in_range(x[2], 0.0, M_PI_2)
             || !in_range(x[3], -M_PI, M_PI)) {
            std::cout << "<<ERROR>> configuration " <<  cart << " -> " << x << " failed" << std::endl;
            ++n;
          }

#if 0
          std::cout << "<<INFO>> checked configuration "
                    << cart << "\r";
#endif
        }
      }
    }
  }

  if(n > 0) {
    std::cerr << "<<ERROR>> " << n << " tests of conversion CARTESIAN -> HALF failed" << std::endl;
  }

  return n;
}
