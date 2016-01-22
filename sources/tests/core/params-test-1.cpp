/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
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

#include <cmath>
#include <cstdlib>
#include <iostream>

/* Test different configurations for the Half / Cartesian parametrization.  */
int main(int argc, char** argv) {

  // Number of failed tests
  unsigned int n = 0, total = 0;

  const int K = 100;
  const int L = 100;
  for(int k=0; k<=K; ++k) {
    for(int l=0; l<=L; ++l) {
      vec x(2), xx(2), y(3);
      vec cart(6);

      x[0] = double(k) / double(K);
      x[1] = double(l) / double(L);

      if(x[0]*x[0] + x[1]*x[1] > 1.0)
        continue;

      total++;
      params::convert(&x[0], params::STARK_2D, params::CARTESIAN, &cart[0]);
      params::convert(&cart[0], params::CARTESIAN, params::RUSIN_TH_TD_PD, &y[0]);
      params::convert(&y[0], params::RUSIN_TH_TD_PD, params::STARK_2D, &xx[0]);

      if(!close_to(x[0], xx[0]) || !close_to(x[1], xx[1])) {
        std::cerr << "x  = " << x << std::endl;
        std::cerr << "c  = " << cart << std::endl;
        std::cerr << "y  = " << y << std::endl;
        std::cerr << "xx = " << xx << std::endl;
        ++n;
      }
    }
  }

  if(n > 0) {
      std::cerr << "<<ERROR>> " << n << " out of " << K * L
                << " conversions tests STARK -> CARTESIAN "
                << "-> HALF -> STARK failed" << std::endl;
  }

  return n > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
