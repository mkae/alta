// ALTA includes
#include <core/params.h>

// STL includes
#include <cmath>
#include <iostream>

bool closeTo(double a, double b) {
	return std::abs(a-b) < 1.0E-10;
}

/* Test different configurations for the Half / Cartesian parametrization
 * Returns: 0 is every test passes
 *          n > 0 when n tests did no pass
 */
int main(int argc, char** argv) {

	// Number of failed tests
	int n = 0;

	// Test the theta_h = 0 configuration
	// This correspond to wi = R(wo), where R is the reflection with respect
	// to the normal vector (0,0,1). 
	for(int k=0; k<=90; ++k) {
		const double theta = 0.5 * M_PI * double(k) / 90.0;
		vec cart(6);
		cart[0] = sin(theta);
		cart[1] = 0.0;
		cart[2] = cos(theta);
		cart[3] =-sin(theta);
		cart[4] = 0.0;
		cart[5] = cos(theta);
		vec x(3);

		params::convert(&cart[0], params::CARTESIAN, params::RUSIN_TH_TD_PD, &x[0]);

		if(!closeTo(x[0], 0.0) || !closeTo(x[1], theta) || !closeTo(x[2], 0.0)) {
			std::cout << "<<ERROR>> configuration " <<  cart << " -> " << x << " failed" << std::endl;
			++n;
		}
	}

	// Test the configuration theta_H = 0 by converting half/diff vector
	// to cartesian. This should convert to wi = R(wo) configuration.
	for(int k=0; k<=90; ++k) {
		const double theta = 0.5 * M_PI * double(k) / 90.0;
		vec half(6);
		half[0] = 0.0;
		half[1] = theta;
		half[2] = 0.0;
		vec cart(6);

		params::convert(&half[0], params::RUSIN_TH_TD_PD, params::CARTESIAN, &cart[0]);

		if(!closeTo(cart[2], cos(theta)) || !closeTo(cart[2], cos(theta))) {
			std::cout << "<<ERROR>> configuration " <<  half << " -> " << cart << " failed" << std::endl;
			++n;
		}
	}

	if(n > 0) {
		std::cerr << "<<ERROR>> " << n << " tests of conversion CARTESIAN -> HALF failed" << std::endl;
	}

	return n;
}