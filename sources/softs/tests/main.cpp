#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QPluginLoader>
#include <QtPlugin>
#include <QCoreApplication>
#include <QDir>
#include <QTime>

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
	QCoreApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	plugins_manager manager(args) ;

	int nb_tests_failed = 0;

	// Parametrization tests
	//
	nb_tests_failed += parametrization_tests();

	// Evaluation tests
	//


	std::cout << "<<INFO>> " << nb_tests_failed << " tests failed" << std::endl;


	return nb_tests_failed;
}

int parametrization_tests()
{
	// Params
	int nb_tests_failed = 0;


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
	std::cout << "<<DEBUG>> rusin 3d: " << rusi << std::endl;
	std::cout << "<<DEBUG>> cartesian: " << spherical << std::endl;
	std::cout << std::endl;


	// Convert issue #1
	vec cart2(6);
	spherical[0] = 0; spherical[1] = 0; spherical[2] = 1.51844; spherical[3] = -2.96706;
	params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart2[0]);
	std::cout << "<<DEBUG>> spherical before: " << spherical << std::endl;
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
	std::cout << "<<DEBUG>> rusin after: " << rusi << std::endl;
	std::cout << "<<DEBUG>> spherical after: " << spherical << std::endl;
	std::cout << "<<DEBUG>> " << cart << " / " << cart2 << std::endl;
	{
	double dot  = cart[0]*cart[3] + cart[1]*cart[4] + cart[2]*cart[5];
	double dot2 = cart2[0]*cart2[3] + cart2[1]*cart2[4] + cart2[2]*cart2[5];
	std::cout << "<<DEBUG>> dot = " << dot << ", " << "dot2 = " << dot2 << std::endl;
	}
	std::cout << std::endl;

	return nb_tests_failed;
}
