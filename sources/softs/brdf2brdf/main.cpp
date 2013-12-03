/*! \package brdf2brdf
 *  \ingroup commands
 *  \brief
 *  This command exports a \ref function object to software specific file.
 *  \details
 *  <h3>Parameters</h3>
 *  <ul>
 *		<li><b>\-\-input <i>filename</i></b> ALTA function file to be loaded.
 *		</li>
 *  </ul>
 */
#include <core/args.h>
#include <core/function.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;
	
	if(args.is_defined("help")) {
		std::cout << "Usage: brdf2brdf [options] --input data.file --output data.file" << std::endl ;
		std::cout << "Re-export a function object to another output format."<< std::endl ;
		std::cout << std::endl;
		std::cout << "Mandatory arguments:" << std::endl;
		std::cout << "  --input    [filename]" << std::endl;
		std::cout << "  --output   [filename]" << std::endl;
		std::cout << std::endl;
		std::cout << "Optional arguments:" << std::endl;
		std::cout << "  --export   [type]      Name of the export format used to save the outputed" << std::endl ;
		std::cout << "                         function file. Available types are: alta, matlab," << std::endl ;
		std::cout << "                         explorer or shader." << std::endl ;
		return 0 ;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}

    // Load the function
    function* f = plugins_manager::get_function(args["input"]);

    // Save it
    f->save(args["output"], args) ;

	return 0 ;
}
