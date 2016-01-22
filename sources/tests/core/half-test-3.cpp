/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// ALTA includes
#include <core/params.h>
#include <tests.h>

using namespace alta;
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
      auto fail = [&](const vec& erroneous_result)
      {
          std::cerr << "FAIL: conversion to Rusin. "
                    << erroneous_result.size() << "D: "
                    << "θ = " << theta_i
                    << " φ = " << phi_i
                    << " → " << erroneous_result << std::endl;
          n++;
      };

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
          fail(x);

      vec y(4);

      // Likewise for the 3D parametrization.
      params::convert(&cart[0], params::CARTESIAN,
                      params::RUSIN_TH_TD_PD, &y[0]);
      if (!close_to(y[0], theta_h)
          || !close_to(y[1], theta_d)
          || !close_to(y[2], phi_d))
          fail(y);

      // And 2D.
      params::convert(&cart[0], params::CARTESIAN,
                      params::RUSIN_TH_TD, &y[0]);
      if (!close_to(y[0], theta_i)
          || !close_to(y[1], 0.))
          fail(y);

#if 0                                       // XXX: currently not implemented
      vec z(6);
      params::convert(&cart[0], params::CARTESIAN,
                      params::RUSIN_VH_VD, &z[0]);
      if (!close_to(z[0], cart[0])                // h⃗
          || !close_to(z[1], cart[1])
          || !close_to(z[2], cart[2])
          || !close_to(z[3], 0.)                  // D⃗
          || !close_to(z[4], 0.)
          || !close_to(z[5], 0.))
          fail(z);
#endif

      // Lastly, 1D.
      params::convert(&cart[0], params::CARTESIAN,
                      params::RUSIN_VH, &y[0]);
      if (!close_to(y[0], cart[0])                // h⃗
          || !close_to(y[1], cart[1])
          || !close_to(y[2], cart[2]))
          fail(y);
  }}

  if (n > 0)
      std::cerr << n << " failures" << std::endl;


  return n > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
