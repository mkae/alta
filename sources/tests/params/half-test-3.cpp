/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Inria

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
#include <cassert>
#include <cstdlib>


int main(int argc, char** argv) {

  // Number of failed tests
  unsigned int n = 0;

  const int step = 10;                    // degrees to add at each iteration

  // Check conversion from Cartesian coordinates to Rusinkiewicz parameters
  // when ωi = ωo = h⃗.
  for (auto&& theta_i : angle_range<double>(0, 90, step)) {
  for (auto&& phi_i : angle_range<double>(-180, 180, step)) {

      if (theta_i == 0. && phi_i != 0.)
          // It doesn't make sense for φi to be different from zero here.
          continue;
      if (phi_i == 0. && theta_i != 0.)
          // Likewise.
          continue;

      auto theta_o = theta_i;
      auto phi_o = phi_i;

      vec cart(6);
      cart[0] = cos(phi_i)*sin(theta_i);
      cart[1] = sin(phi_i)*sin(theta_i);
      cart[2] = cos(theta_i);
      cart[3] = cos(phi_o)*sin(theta_o);
      cart[4] = sin(phi_o)*sin(theta_o);
      cart[5] = cos(theta_o);
      vec x(4);

      params::convert(&cart[0], params::CARTESIAN,
                      params::RUSIN_TH_PH_TD_PD, &x[0]);

      auto theta_h = x[0];
      auto phi_h = x[1];
      auto theta_d = x[2];
      auto phi_d = x[3];

#if 0
      std::cout << "conf " << theta_i << " " << phi_i
                << " → φd = " << phi_d
                << std::endl;
#endif

      // Make sure θh = θi and φh = φi, and θd = φd = 0.
      if (!close_to(theta_h, theta_i)
          || !close_to(phi_h, phi_i)
          || !close_to(theta_d, 0.)
          || !close_to(phi_d, 0.))
      {
          std::cerr << "FAIL: θ = " << theta_i
                    << " φ = " << phi_i
                    << " → " << x << std::endl;
          n++;
      }
  }}

  if (n > 0)
      std::cerr << n << " failures" << std::endl;


  return n > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
