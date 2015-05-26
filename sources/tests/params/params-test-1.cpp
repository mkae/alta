// ALTA includes
#include <core/params.h>

// STL includes
#include <cmath>
#include <iostream>

bool closeTo(double a, double b) {
	return std::abs(a-b) < 1.0E-6;
}

/* Test different configurations for the Half / Cartesian parametrization
 * Returns: 0 if every test passes
 *          n > 0 when n tests did no pass
 */
int main(int argc, char** argv) {

	// Number of failed tests
	int n = 0;
	
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

			params::convert(&x[0], params::STARK_2D, params::CARTESIAN, &cart[0]);
			params::convert(&cart[0], params::CARTESIAN, params::RUSIN_TH_TD_PD, &y[0]);
			params::convert(&y[0], params::RUSIN_TH_TD_PD, params::STARK_2D, &xx[0]);

			if(!closeTo(x[0], xx[0]) || !closeTo(x[1], xx[1])) {
				std::cerr << "x  = " << x << std::endl;
				std::cerr << "c  = " << cart << std::endl;
				std::cerr << "y  = " << y << std::endl;
				std::cerr << "xx = " << xx << std::endl;
				++n;
			}
		}
	}

	if(n > 0) {
		std::cerr << "<<ERROR>> " << n << " tests of conversion CARTESIAN -> HALF failed" << std::endl;
	}

	return n;
}
