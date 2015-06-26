/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
   Copyright (C) 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// ALTA includes
#include <core/params.h>

// STL includes
#include <cmath>
#include <iostream>

bool closeTo(double a, double b) {
	return std::abs(a-b) < 1.0E-10;
}

static bool inRange(double x, double a, double b) {
	return x >= a && x <= b;
}

/* Test different configurations for the Half / Cartesian parametrization
 * Returns: 0 is every test passes
 *          n > 0 when n tests did no pass
 */
int main(int argc, char** argv) {

	// Number of failed tests
	int n = 0;

	const double d2r = M_PI / 180.0;
	const int nphi   = 360;
	const int ntheta = 90;
	const int step   = 10;

	// Sample the full hemisphere for both ωi and ωo.  Check if the converted
	// vector is within the range of the parametrization.  θh and θd should be
	// within [0..π/2], and Φd should be within [-π..π].
	for(int ti=0; ti<=ntheta; ti+=step) {
		for(int pi=0; pi<=nphi; pi+=step) {

			const double theta_i = d2r*double(ti);
			const double phi_i   = d2r*double(pi);

			#pragma omp parallel for
			for(int to=0; to<=ntheta; to+=step) {
				for(int po=0; po<=nphi; po+=step) {

					const double theta_o = d2r*double(to);
					const double phi_o   = d2r*double(po);

					vec cart(6);
					cart[0] = cos(phi_i)*sin(theta_i);
					cart[1] = sin(phi_i)*sin(theta_i);
					cart[2] = cos(theta_i);
					cart[3] = cos(phi_o)*sin(theta_o);
					cart[4] = sin(phi_o)*sin(theta_o);
					cart[5] = cos(theta_o);
					vec x(4);

					params::convert(&cart[0], params::CARTESIAN, params::RUSIN_TH_PH_TD_PD, &x[0]);

					#pragma omp critical (n)
					if(!inRange(x[0], 0.0, M_PI_2) || !inRange(x[2], 0.0, M_PI_2)
             || !inRange(x[3], -M_PI, M_PI)) {
						std::cout << "<<ERROR>> configuration " <<  cart << " -> " << x << " failed" << std::endl;
						++n;
					}

          if (po % 10000 == 0) {
              std::cout << "<<INFO>> Check configuration " << po + nphi*(to + ntheta*(pi + nphi*ti)) << " / " << ntheta*ntheta*nphi*nphi << "      \r";
          }
				}
			}
		}
	}

	if(n > 0) {
		std::cerr << "<<ERROR>> " << n << " tests of conversion CARTESIAN -> HALF failed" << std::endl;
	}

	return n;
}
