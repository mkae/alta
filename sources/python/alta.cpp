/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// Boost includes
#include <boost/python.hpp>

// ALTA include
#include <core/common.h>
#include <core/ptr.h>
#include <core/function.h>
#include <core/plugins_manager.h>

// STL include
#include <iostream>

#define bp boost::python


/* The following code register ALTA's shared pointer as a valid shared ptr
 * to be used by boost::python .
 */
template <typename T> 
T* get_pointer(ptr<T> const& p) {
  return const_cast<T*>(p.get());
}

namespace boost {
	namespace python {
   		template <typename T>
    	struct pointee< ::ptr<T> > {
        	typedef T type;
    	};
	}
}


/* Wrapper to ALTA's vec class. This is only here to allow init with Python's
 * list.
 *
 * TODO: Make sure that the value passed to this vector are floatting point
 *       convertible.
 */
struct python_vec : public vec {
    python_vec() : vec() {}
	python_vec(const vec& x) : vec(x) {}
    python_vec(const bp::list& l) : vec(bp::len(l)) {
        for(auto i=0; i<bp::len(l); ++i) {
            (*this)[i] = bp::extract<double>(l[i]);
        }
    }
};
std::ostream &operator<<(std::ostream &out, const python_vec &x) {
	out << "[";
	for(int i=0; i<x.size(); ++i) {
		if(i != 0) { out << ", "; }
		out << x[i]; 
	}
	return out << "]";
}



/* This class is a wrapper to ALTA's arguments class to add Python specific
 * behaviour such as dictionnary initialization.
 *
 * Right now it does not handle automatic conversion in function call for 
 * example. The following code is not possible:
 *
 *    import alta
 *    alta.get_data('data_merl', {'params' : 'STARK_2D'})
 *
 * Instead, one has to construct the arguments object from the ALTA library
 * to use it afterwards:
 *
 *    import alta
 *    args = alta.arguments({'params' : 'STARK_2D'})
 *    alta.get_data('data_merl', args)
 */
struct python_arguments : public arguments {
	python_arguments() : arguments() {}
	python_arguments(bp::dict d) : arguments() {
		bp::list keys = d.keys();
		for(int i=0; i<bp::len(keys); ++i) {
			const std::string s_key = bp::extract<std::string>(keys[i]);
			const std::string s_val = bp::extract<std::string>(d[keys[i]]);
			this->update(s_key, s_val);
		}
	}
};

/* Create a data object from a plugin's name and the data filename. This 
 * function is here to accelerate the loading of data file.
 */
ptr<data> load_data(const std::string& plugin_name, const std::string& filename) {
	ptr<data> d = plugins_manager::get_data(plugin_name);
	d->load(filename);
	return d;
}
ptr<data> get_data_with_args(const std::string& plugin_name, const python_arguments& args) {
	return plugins_manager::get_data(plugin_name, args);
}
ptr<data> get_data(const std::string& plugin_name) {
	return plugins_manager::get_data(plugin_name);
}

/* Creating functions for the plugins_manager calls
 * 
 * TODO: Throw python exceptions if the function is not correctly created.
 *       Those function should disapear when the return type of get_Function
 *       in the plugin_manager will be ptr<function>.
 */
ptr<function> get_function(const std::string& plugin_name) {
    ptr<function> func(plugins_manager::get_function(plugin_name));
    if(!func) {
    	std::cerr << "<<ERROR>> no function created" << std::endl;
    }
	return func;
}
ptr<function> get_function_from_args(const python_arguments& args) {
    ptr<function> func(plugins_manager::get_function(args));
    if(!func) {
    	std::cerr << "<<ERROR>> no function created" << std::endl;
    }
    return func;
}


/* Softs functions. Those function recopy the softs main function, without
 * the command line arguments.
 * TODO: Add the command line arguments in the parameters
 */
void data2data(const data* d_in, data* d_out) {
	#pragma omp parallel for
	for(int i=0; i<d_out->size(); ++i)
	{
		vec temp(d_in->dimX());
		vec cart(6);
		vec y(d_in->dimY());

		// Copy the input vector
		vec x = d_out->get(i);
		params::convert(&x[0], d_out->parametrization(), params::CARTESIAN, &cart[0]);

		if(cart[2] >= 0.0 || cart[5] >= 0.0) {
			params::convert(&cart[0], params::CARTESIAN, d_in->parametrization(), &temp[0]);
			y = d_in->value(temp);
		} else {
			y.setZero();
		}

		params::convert(&y[0], d_in->output_parametrization(), d_in->dimY(), d_out->output_parametrization(), d_out->dimY(), &x[d_out->dimX()]);

		d_out->set(x);
	}	
}

void fit_data(ptr<fitter>& f, const ptr<data>& d, ptr<function>& fn, const arguments& args) {
	f->fit_data(d, fn, args);
}


/* Exporting the ALTA module 
 */
BOOST_PYTHON_MODULE(alta)
{
	// Argument class
	//
	bp::class_<python_arguments>("arguments")
		.def(bp::init<>())
		.def(bp::init<bp::dict>())
		.def("update", &arguments::update);


	// Vector class
	//
	// TODO: There is a conversion issue right now that prevents us from using vectors
	// within Python. This needs to be investiguated.
	bp::class_<python_vec>("vec")
		.def(bp::init<vec>())
		.def(bp::init<bp::list>())
		.def(bp::self_ns::str(bp::self_ns::self))
		/*.def(" __setitem__", &my_vec::operator[])*/;
	bp::implicitly_convertible<vec, python_vec>();
	bp::implicitly_convertible<python_vec, vec>();


	// Function interface
	//
	bp::class_<function, ptr<function>, boost::noncopyable>("function", bp::no_init)
		.def("value", &function::value)
		.def("load", &function::load)
		.def("save",  &function::save);
	bp::def("get_function", get_function);
	bp::def("get_function", get_function_from_args);



	// Data interface
	//
	bp::class_<data, ptr<data>, boost::noncopyable>("data", bp::no_init)
		.def("size", &data::size)
		.def("get",  &data::get)
		//.def("set",  &data::set)
		.def("load", static_cast< void(data::*)(const std::string&)>(&data::load))
		.def("save", &data::save);
	bp::def("get_data",  get_data);
	bp::def("get_data",  get_data_with_args);
	bp::def("load_data", load_data);


	// Fitter interface
	//
	bp::class_<fitter, ptr<fitter>, boost::noncopyable>("fitter", bp::no_init)
		.def("fit_data", &fitter::fit_data);
	bp::def("get_fitter", plugins_manager::get_fitter);
	bp::def("fit_data",   fit_data);

	// Softs
	//
	bp::def("data2data", data2data);	
}
