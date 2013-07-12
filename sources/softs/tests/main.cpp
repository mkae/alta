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

	
	// Test Rusinkevich parametrization
	vec cart(6);
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


	return nb_tests_failed;
}
