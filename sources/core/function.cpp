#include "function.h"		

#include "common.h"
#include "plugins_manager.h"

void function::bootstrap(const data*, const arguments& args)
{
    // If the bootstrap option contains a filename, load it
    if(args.is_defined("bootstrap") && !args["bootstrap"].empty())
    {
        std::ifstream in(args["bootstrap"].c_str());
        if(in.is_open())
        {
            load(in);
        }
    }
}

// Acces to the domain of definition of the function
void function::setMin(const vec& min)
{
#ifdef DEBUG
    assert(min.size() == _nX) ;
#endif
    _min = min ;
}
void function::setMax(const vec& max)
{
#ifdef DEBUG
    assert(max.size() == _nX) ;
#endif
    _max = max ;
}
vec function::getMin() const
{
    return _min ;
}
vec function::getMax() const
{
    return _max ;
}

void function::save(const std::string& filename, const arguments& args) const
{
    bool is_cpp      = args["export"] == "C++";
    bool is_explorer = args["export"] == "explorer";
    bool is_shader   = args["export"] == "shader" || is_explorer;
    bool is_matlab   = args["export"] == "matlab";

	// Open the file
	std::ofstream file(filename.c_str());
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open output file for writing" << std::endl;
	}

    if(is_explorer)
    {
        file << "analytic" << std::endl;
        file << std::endl;
        file << "::begin shader" << std::endl;
    }

	// Save common header
	save_header(file, args);
    file << std::endl;

	// Save function definition
	save_body(file, args);
    file << std::endl;

    if(is_cpp)
    {
        file << "vec brdf(const vec& in, const vec& file)" << std::endl;
        file << "{" << std::endl;
        file << "\tvec res(" << dimY() << ");" << std::endl;
        file << "\t";
    }
    else if(is_matlab)
    {
        file << "function res = brdf(in, file)" << std::endl;
        file << "{" << std::endl;
        file << "\tres = zeros(" << dimY() << ");" << std::endl;
        file << "\t";
    }
    else if(is_shader)
    {
        file << "vec3 BRDF(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y)" << std::endl;
        file << "{" << std::endl;
        file << "\tvec3 res = ";
    }

	// Save fit data
    save_call(file, args);
    if(is_cpp || is_shader)
    {
        file << ";" << std::endl;
        file << "\treturn res;" << std::endl;
        file << "}" << std::endl;
    }
    else if(is_matlab)
    {
        file << ";" << std::endl;
        file << "\treturn res;" << std::endl;
        file << "endfunction" << std::endl;
    }
    file << std::endl;

    if(is_explorer)
    {
        file << "::end shader" << std::endl;
    }
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
        //out << "#PARAM_OUT " << params::get_name(output_parametrization()) << std::endl;*
        if(args.is_defined("export-append")) {
            out << args["export-append"] << std::endl;
        }
		out << "#ALTA HEADER END" << std::endl;
		out << std::endl;
	}
}

//! \brief save function specific data. This has no use for ALTA export
//! but allows to factorize the code in the C++ or matlab export by
//! defining function calls that are common to all the plugins.
void function::save_body(std::ostream&, const arguments&) const
{


}

//! \brief save object specific information. For an ALTA export the
//! coefficients will be exported. For a C++ or matlab export, the call
//! to the associated function will be done.
void function::save_call(std::ostream&, const arguments&) const
{

}

//! \brief L2 norm to data.
double function::L2_distance(const data* d) const
{
	double l2_dist = 0.0;
	for(int i=0; i<d->size(); ++i)
	{
		vec dat = d->get(i);
        vec x(d->dimX()), y(d->dimY());
        for(int j=0; j<d->dimX(); ++j) { x[j] = dat[j]; }
        for(int j=0; j<d->dimY(); ++j) { y[j] = dat[d->dimX()+j]; }

		//linf_dist = std::max<double>(linf_dist, std::abs<double>(norm(y-rj->value(dat))));
        l2_dist  += std::pow(norm(y-value(x)), 2);
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
        vec x(d->dimX()), y(d->dimY());
        for(int j=0; j<d->dimX(); ++j) { x[j] = dat[j]; }
        for(int j=0; j<d->dimY(); ++j) { y[j] = dat[d->dimX()+j]; }

		linf_dist = std::max<double>(linf_dist, std::abs(norm(y-value(dat))));
	}

	return linf_dist;
}


void nonlinear_function::bootstrap(const data* d, const arguments& args)
{
    if(args.is_vec("bootstrap"))
    {
        vec p = args.get_vec("bootstrap", nbParameters());
        setParameters(p);
    }
    else
    {
        function::bootstrap(d, args);
    }
}
		
//! \brief get the maximum value for all the parameters in a vector
//! format. The maximum vector and the parameter vector have the same
//! indexing.
vec nonlinear_function::getParametersMax() const
{
	vec M(nbParameters());
	for(int i=0; i<nbParameters(); ++i) 
	{
		M[i] = std::numeric_limits<double>::max();
	}
	return M;
}

//! \brief get the minimum value for all the parameters in a vector
//! format. The minimum vector and the parameter vector have the same
//! indexing.
vec nonlinear_function::getParametersMin() const
{
	vec m(nbParameters());
	for(int i=0; i<nbParameters(); ++i) 
	{
		m[i] = -std::numeric_limits<double>::max();
	}
	return m;
}

void compound_function::bootstrap(const ::data* d, const arguments& args)
{
    const bool global_bootstrap = args.is_defined("bootstrap");

    // Check if the bootstrap is global
    if(global_bootstrap)
    {
        if(args.is_vec("bootstrap"))
        {
            vec p = args.get_vec("bootstrap", nbParameters());
            setParameters(p);
        }
        else
        {
            std::ifstream file;
            file.open(args["bootstrap"].c_str()) ;
            if(file.is_open())
            {
                // Skip the header
                std::string line ;
                do
                {
                    std::getline(file, line) ;
                }
                while(line != "#ALTA HEADER END");

                // Load the file
                this->load(file);
            }
            else
            {
                std::cerr << "<<ERROR>> you must provide a vector of parameters or a file to load with the bootstrap" << std::endl;
            }
        }
    }

    // Per function bootstrap
    for(unsigned int i=0; i<fs.size(); ++i)
    {
        if(!global_bootstrap)
        {
            fs[i]->bootstrap(d, fs_args[i]);
        }
    }
}
