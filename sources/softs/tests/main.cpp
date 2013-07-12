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

#define EPSILON 1.0E-10

int main(int argc, char** argv)
{
	QCoreApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	plugins_manager manager(args) ;

	int nb_tests_failed = 0;

	// Parametrization tests
	//
	

	// Evaluation tests
	//


	std::cout << "<<INFO>> " << nb_tests_failed << " tests failed" << std::endl;


	return nb_tests_failed;
}
