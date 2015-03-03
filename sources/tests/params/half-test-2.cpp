// ALTA includes
#include <core/params.h>

// STL includes
#include <cmath>
#include <iostream>

bool closeTo(double a, double b) {
	return std::abs(a-b) < 1.0E-10;
}

bool inRange(double x, double a, double b) {
	return x >= a && x <= b;
}

/* Test different configurations for the Half / Cartesian parametrization
 * Returns: 0 is every test passes
 *          n > 0 when n tests did no pass
 */
int main(int argc, char** argv) {

	// Number of failed tests
	int n = 0;

	const double PI2 = 0.5*M_PI;
	const double d2r = M_PI / 180.0;
	const int nphi   = 360;
	const int ntheta = 90;
	const int step   = 10;

	// Sample the full hemisphere for both wi and wo. Check if the converted
	// vector is within the range of the parametrization. Theta_h and Theta_d
	// should be within [0..pi/2] and Phi_d should be within [0..2Pi].
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
					if(!inRange(x[0], 0.0, PI2) || !inRange(x[2], 0.0, PI2)) {
						std::cout << "<<ERROR>> configuration " <<  cart << " -> " << x << " failed" << std::endl;
						++n;
					}

					std::cout << "<<INFO>> Check configuration " << po + nphi*(to + ntheta*(pi + nphi*ti)) << " / " << ntheta*ntheta*nphi*nphi << "      \r";
				}
			}
		}
	}

	if(n > 0) {
		std::cerr << "<<ERROR>> " << n << " tests of conversion CARTESIAN -> HALF failed" << std::endl;
	}

	return n;
}