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


// Test the conversion from Cartesian coordinate to Rusinkiewicz
// parametrizations in the special case of grazing angles, specifically where
// ωi = (1,0,0) and ωo = (0,1,0).
int main(int argc, char** argv) {

  // Number of failed tests
  unsigned int n = 0;

  // Special case.
  vec cart(6);
  cart[0] = 1, cart[1] = 0, cart[2] = 0;          // ωi
  cart[3] = 0, cart[4] = 1, cart[5] = 0;          // ωo

  // What we expect.
  double theta_h = M_PI_2;
  double phi_h = M_PI_4;
  double theta_d = M_PI_4;
  double phi_d = -M_PI_2;
  vec h(3);                                       // h⃗
  h[0] = M_SQRT2 / 2, h[1] = M_SQRT2 / 2, h[2] = 0.;

  auto fail = [&](params::input target, const vec& erroneous_result)
  {
      std::cerr << "FAIL: conversion to " << params::get_name(target) << ": "
                << " ωi|ωo = " << cart
                << " → " << erroneous_result << std::endl;
      n++;
  };

  vec x(4);
  params::convert(&cart[0], params::CARTESIAN,
                  params::RUSIN_TH_PH_TD_PD, &x[0]);

  if (!close_to(x[0], theta_h)
      || !close_to(x[1], phi_h)
      || !close_to(x[2], theta_d)
      || !close_to(x[3], phi_d))
      fail(params::RUSIN_TH_PH_TD_PD, x);

  params::convert(&cart[0], params::CARTESIAN,
                  params::RUSIN_TH_TD_PD, &x[0]);

  if (!close_to(x[0], theta_h)
      || !close_to(x[1], theta_d)
      || !close_to(x[2], phi_d))
      fail(params::RUSIN_TH_TD_PD, x);

  params::convert(&cart[0], params::CARTESIAN,
                  params::RUSIN_TH_TD, &x[0]);

  if (!close_to(x[0], theta_h)
      || !close_to(x[1], theta_d))
      fail(params::RUSIN_TH_TD, x);

  params::convert(&cart[0], params::CARTESIAN,
                  params::RUSIN_VH, &x[0]);
  if (!close_to(x[0], h[0])
      || !close_to(x[1], h[1])
      || !close_to(x[2], h[2]))
      fail(params::RUSIN_VH, x);

  if (n > 0)
      std::cerr << n << " failures" << std::endl;


  return n > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
