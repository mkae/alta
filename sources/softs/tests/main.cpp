
/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/param.h>

int main(int argc, char** argv) {

	spherical sphericalParam;
	isotropic_tl_tv isotropicParam;

	vec in(4);
	in[0] = 0.0;
	in[1] = 0.0;
	in[2] = 0.0;
	in[3] = 0.0;

	vec out(2);
	sphericalParam.convert_to(in, isotropicParam, out);

	return 0;
}


#ifdef OLD


#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

#define EPSILON 1.0E-5

bool is_close(double x, double y)
{
	return std::abs(x - y) < EPSILON;
}

int parametrization_tests();

int main(int argc, char** argv)
{

	int nb_tests_failed = 0;

	// Parametrization tests
	//
	nb_tests_failed += parametrization_tests();

	// Evaluation tests
	//

	std::cout << "<<INFO>> " << nb_tests_failed << " tests failed" << std::endl;


	return nb_tests_failed;
}

bool is_close(const vec& a, const vec& b)
{
	double dist = 0.0;
	for(int i=0; i<a.size(); ++i)
	{
		dist += abs(a[i] - b[i]);
	}
	return dist < 1.0E-10;
}

int parametrization_tests()
{
	// Params
	int nb_tests_failed = 0;

	
	// Test RUSIN_TH_TD
	vec rhd(2), cart(6), res(6), vh(3);

	// Convert RUSIN_TH_TD (0,0) to CARTESIAN (0,0,1,0,0,1)
	rhd[0] = 0;	rhd[1] = 0;
	res[0] = 0; res[1] = 0; res[2] = 1; res[3] = 0; res[4] = 0; res[5] = 1;
	params::convert(&rhd[0], params::RUSIN_TH_TD, params::CARTESIAN, &cart[0]);
	if(!is_close(cart, res)) {
		nb_tests_failed++;
		std::cout << "<<FAILED>> out = " << cart << ", while attending " << res << std::endl;
	}

	// Convert RUSIN_TH_TD (0,pi/2) to CARTESIAN (1,0,0,-1,0,0)
	rhd[0] = 0;	rhd[1] = M_PI*0.5;
	res[0] = 1; res[1] = 0; res[2] = 0; res[3] = -1; res[4] = 0; res[5] = 0;
	params::convert(&rhd[0], params::RUSIN_TH_TD, params::CARTESIAN, &cart[0]);
	if(!is_close(cart, res)) {
	  	nb_tests_failed++;
		std::cout << "<<FAILED>> out = " << cart << ", while attending " << res << std::endl;
	}

	// Convert RUSIN_TH_TD (0, pi/2) to CARTESIAN (1,0,0,-1,0,0) and to RUSIN_VH (0,0,1)
	rhd[0] = 0;	rhd[1] = M_PI*0.5;
	res[0] = 1; res[1] = 0; res[2] = 0; res[3] = -1; res[4] = 0; res[5] = 0;
	params::convert(&rhd[0], params::RUSIN_TH_TD, params::CARTESIAN, &cart[0]);
	if(!is_close(cart, res)) {
	  	nb_tests_failed++;
		std::cout << "<<FAILED>> out = " << cart << ", while attending " << res << std::endl;
	}
	else
	{
		params::convert(&rhd[0], params::RUSIN_TH_TD, params::RUSIN_VH, &vh[0]);
		res.resize(3); res[0] = 0; res[1] = 0; res[2] = 1;
		if(!is_close(vh, res)) {
			nb_tests_failed++;
			std::cout << "<<FAILED>> out = " << vh << ", while attending " << res << std::endl;
		}
	}
/*

	// Test the rotation code
	vec x(3);
	x[0] = 0; x[1] = 0; x[2] = 1;
	params::rotate_binormal(&x[0], 0.5*M_PI);
	std::cout << "<<DEBUG>> x = " << x << std::endl;

	params::rotate_binormal(&x[0], -0.5*M_PI);
	std::cout << "<<DEBUG>> x = " << x << std::endl;
	std::cout << "<<DEBUG>> acos(x[2]) = " << acos(x[2]) << std::endl;
	
	// Test Rusinkevich parametrization
	vec cart(6), spherical(4);
	vec rusi(3);

	// Equal directions test, when the PHI_D is ZERO
	rusi[0] = 0.25*M_PI; rusi[1] = 0.0; rusi[2] = 0.0;
	params::convert(&rusi[0], params::RUSIN_TH_TD_PD, params::CARTESIAN, &cart[0]);

	const double dot = cart[0]*cart[3] + cart[1]*cart[4] + cart[2]*cart[5];

	if(!is_close(cart[0], cart[3]) || !is_close(cart[1], cart[4]) ||
	   !is_close(cart[2], cart[5]) || !is_close(dot, 1.0))
	{
		std::cout << "<<DEBUG>> rusin 3d: " << rusi << std::endl;
		std::cout << "<<DEBUG>> cartesian: " << cart << std::endl;
		std::cout << "<<DEBUG>> dot: " << dot << std::endl;
		std::cout << std::endl;
		nb_tests_failed++;
	}


	// Pathological case when THETA_H and THETA_D are equal to
	// PI/4 and PHI_D PI, the conversion seems to fail.
	rusi[0] = 0.759218; rusi[1] = 0.759218; rusi[2] = 3.14159;
	try
	{
		params::convert(&rusi[0], params::RUSIN_TH_TD_PD, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
		params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::RUSIN_TH_TD_PD, &rusi[0]);
	}
	catch(...)
	{
		std::cout << "<<ERROR>> the conversion failed" << std::endl;
		std::cout << "<<DEBUG>> rusin 3d: " << rusi << std::endl;
		std::cout << "<<DEBUG>> cartesian: " << spherical << std::endl;
		nb_tests_failed++;
    }


	// Convert issue #1
	vec cart2(6);
	spherical[0] = 0; spherical[1] = 0; spherical[2] = 1.51844; spherical[3] = -2.96706;
	params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart2[0]);
	try
	{
		params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::RUSIN_TH_TD_PD, &rusi[0]);
		params::convert(&rusi[0], params::RUSIN_TH_TD_PD, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
		params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart[0]);
	}
	catch(...)
	{
		std::cout << "<<ERROR>> the conversion failed" << std::endl;
		std::cout << "<<DEBUG>> rusin 3d: " << rusi << std::endl;
		std::cout << "<<DEBUG>> spherical: " << spherical << std::endl;
		nb_tests_failed++;
    }



	/// Test on a known couple of directions
	/// in = [0, 0, 1] out = [1, 0, 0]

	// Convert from Cartesian to spherical
	cart[0] = 0; cart[1] = 0; cart[2] = 1; 
	cart[3] = 1; cart[4] = 0; cart[5] = 0; 
    params::convert(&cart[0], params::CARTESIAN, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
	 std::cout << "From cartesian to spherical" << std::endl;
	 std::cout << spherical << std::endl << std::endl;
    
	 params::convert(&cart[0], params::CARTESIAN, params::RUSIN_TH_TD_PD, &rusi[0]);
	 std::cout << "From cartesian to rusi" << std::endl;
    std::cout << rusi << std::endl << std::endl;
    
	 params::convert(&rusi[0], params::RUSIN_TH_TD_PD, params::CARTESIAN, &cart[0]);
	 std::cout << "From rusi to cartesian" << std::endl;
    std::cout << cart << std::endl << std::endl;

    params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::RUSIN_TH_TD_PD, &rusi[0]);
	 std::cout << "From spherical to rusi" << std::endl;
    std::cout << rusi << std::endl << std::endl;

    params::convert(&rusi[0], params::RUSIN_TH_TD_PD, params::SPHERICAL_TL_PL_TV_PV, &spherical[0]);
	 std::cout << "From rusi to spherical" << std::endl;
    std::cout << spherical << std::endl << std::endl;

    params::convert(&rusi[0], params::RUSIN_TH_TD_PD, params::CARTESIAN, &cart[0]);
	 std::cout << "From rusi to cartesian" << std::endl;
    std::cout << cart << std::endl << std::endl;

    params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::RUSIN_TH_TD_PD, &rusi[0]);
	 std::cout << "From spherical to rusi" << std::endl;
    std::cout << rusi << std::endl << std::endl;
*/
	return nb_tests_failed;
}

#endif
