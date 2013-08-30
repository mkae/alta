#include "function.h"		

#include "common.h"

void function::save(const std::string& filename, const arguments& args) const
{
	// Open the file
	std::ofstream file(filename);
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open output file for writing" << std::endl;
	}

	// Save common header
	save_header(file, args);

	// Save function definition
	save_body(file, args);

	// Save fit data
	save_call(file, args);
}
		
//! \brief save the header of the output function file. The header should
//! store general information about the fit such as the command line used
//! the dimension of the fit. L2 and L_inf distance could be added here.
void function::save_header(std::ostream& out, const arguments& args) const
{
	if(!args.is_defined("export"))
	{
		out << "#ALTA FUNC HEADER" << std::endl;
		out << "#CMD " << args.get_cmd() << std::endl;
		out << "#DIM " << _nX << " " << _nY << std::endl;
		out << "#PARAM_IN  " << params::get_name(input_parametrization()) << std::endl;
		//out << "#PARAM_OUT " << params::get_name(output_parametrization()) << std::endl;
		out << "#ALTA HEADER END" << std::endl;
		out << std::endl;
	}
}

//! \brief save function specific data. This has no use for ALTA export
//! but allows to factorize the code in the C++ or matlab export by
//! defining function calls that are common to all the plugins.
void function::save_body(std::ostream& out, const arguments& args) const
{
	bool is_cpp    = args["export"] == "C++";
	bool is_shader = args["export"] == "shader";
	bool is_matlab = args["export"] == "matlab";

	if(is_cpp)
	{
		out << "vec brdf(const vec& in, const vec& out)" << std::endl;
		out << "{" << std::endl;
		out << "\tvec res(" << dimY() << ");" << std::endl;
	}
	else if(is_matlab)
	{
		out << "function res = brdf(in, out)" << std::endl;
		out << "\tres = zeros(" << dimY() << ");" << std::endl;
	}
	else if(is_shader)
	{
		out << "vec3 brdf(vec3 in, vec3 out)" << std::endl;
		out << "\tvec3 res = vec3(0.0f);" << std::endl;
	}
}

//! \brief save object specific information. For an ALTA export the
//! coefficients will be exported. For a C++ or matlab export, the call
//! to the associated function will be done.
void function::save_call(std::ostream& out, const arguments& args) const
{
	bool is_cpp    = args["export"] == "C++";
	bool is_shader = args["export"] == "shader";
	bool is_matlab = args["export"] == "matlab";

	if(is_cpp || is_shader)
	{
		out << "\treturn res;" << std::endl;
		out << "}" << std::endl;
	}
	else if(is_matlab)
	{
		out << "endfunction" << std::endl;
	}
}

//! \brief L2 norm to data.
double function::L2_distance(const data* d) const
{
	double l2_dist = 0.0;
	for(int i=0; i<d->size(); ++i)
	{
		vec dat = d->get(i);
		vec y(d->dimY());
		for(int j=0; j<d->dimY(); ++j)
			y[j] = dat[d->dimX()+j];

		//linf_dist = std::max<double>(linf_dist, std::abs<double>(norm(y-rj->value(dat))));
		l2_dist  += std::pow(norm(y-value(dat)), 2);
	}
	l2_dist = std::sqrt(l2_dist / d->size());
	return l2_dist;
}

//! \brief Linf norm to data.
double function::Linf_distance(const data* d) const
{

	double linf_dist = 0.0;
	for(int i=0; i<d->size(); ++i)
	{
		vec dat = d->get(i);
		vec y(d->dimY());
		for(int j=0; j<d->dimY(); ++j)
			y[j] = dat[d->dimX()+j];

		linf_dist = std::max<double>(linf_dist, std::abs(norm(y-value(dat))));
	}

	return linf_dist;
}
