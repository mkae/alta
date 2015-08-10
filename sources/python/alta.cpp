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
#include <core/rational_function.h>
#include <core/plugins_manager.h>
#include <core/vertical_segment.h>

// STL include
#include <iostream>

#define bp boost::python


/* The following code register ALTA's shared pointer as a valid shared ptr
 * to be used by boost::python .
 */
template <typename T> T* get_pointer(ptr<T>& p) {
   return const_cast<T*>(p.get());
}

template <typename T> const T* get_pointer(const ptr<T>& p) {
   return p.get();
}

namespace boost {
   namespace python {

      template <typename T> struct pointee< ::ptr<T> > {
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

/* Create special converter from and to list/vec */
struct list_to_vec_converter {

	// Constructor add the converter to the registry
	list_to_vec_converter() {
		bp::converter::registry::push_back(
				&convertible,
				&construct,
				bp::type_id<vec>());
	}

	// Is the Python object convertible to a vec ?
	static void* convertible(PyObject* obj_ptr) {
		if (!PyList_Check(obj_ptr)) return 0;
		return obj_ptr;
	}

	// From a PyObject, construct a vector object
	static void construct(
			PyObject* obj_ptr,
			bp::converter::rvalue_from_python_stage1_data* data) {

		auto size = PyList_Size(obj_ptr);
		vec* _vec = new vec(size);
		for(auto i=0; i<size; ++i) {
			auto pyitem = PyList_GetItem(obj_ptr, i);
			(*_vec)[i] = PyFloat_AsDouble(pyitem);
		}

		data->convertible = (void*) _vec;;
	}
};
struct vec_to_list_converter {

	// From a vec, create a PyObject
	static PyObject* convert(const vec& x) {
		auto obj = PyList_New(x.size());
		for(auto i=0; i<x.size(); ++i) {
			PyList_SetItem(obj, i, PyFloat_FromDouble(x[i]));	
		}

		return obj;
	}
};

/* Specific function call to acces a vector's element
 */
double vec_get_item(const vec& x, int i) {
	return x(i);
}

/* Specific function call to set a vector's element
 */
void vec_set_item(vec& x, int i, double a) {
	x(i) = a;
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

/* Load a function object from a file. The arguments object is never used here.
 * The file is supposed to load the function correctly.
 */
ptr<function> load_function(const std::string& filename) {
   return plugins_manager::load_function(filename);
}
ptr<function> load_function_with_args(const std::string& filename, const arguments&) {
   return plugins_manager::load_function(filename);
}

/* Loading a function object from file
 *
 * TODO: Add exceptions
 */
void load_from_file_with_args(const ptr<function>& func, const std::string& filename,
                              const arguments& args) {

   // Open a stream
   std::ifstream file;
   file.open(filename.c_str()) ;
   if(!file.is_open()) {
      std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
      return;
   }

   // Parse the associated header
   arguments header = arguments::parse_header(file);

   // Load the function stream
   func->load(file);
}
void load_from_file(const ptr<function>& func, const std::string& filename) {

   arguments args;
   load_from_file_with_args(func, filename, args);
}

/* Operators on function object. This provide the ability to create compounds
 * and product in the command line. This is only possible for nonlinear_functions
 * 
 * TODO: The compound and product function should store the shared pointer to the
 * function objects. They might stay in memory longer than the input functions.
 */
ptr<function> add_function(const ptr<function>& f1, const ptr<function>& f2) {
	// Convert to a nonlinear function
	ptr<nonlinear_function> nf1 = dynamic_pointer_cast<nonlinear_function>(f1);
	ptr<nonlinear_function> nf2 = dynamic_pointer_cast<nonlinear_function>(f2);

	arguments args;
	compound_function* cf = new compound_function();

	if(nf1 && nf2) {
		cf->push_back(nf1, args);
		cf->push_back(nf2, args);
	
		return ptr<function>(cf);

	// Failure case, one of the function is a NULL ptr.
	} else {
		std::cerr << "<<ERROR>> One of the input functions is NULL" << std::endl;
		return ptr<function>(NULL);
	}
}
ptr<function> mult_function(const ptr<function>& f1, const ptr<function>& f2) {
	// Convert to a nonlinear function
	ptr<nonlinear_function> nf1 = dynamic_pointer_cast<nonlinear_function>(f1);
	ptr<nonlinear_function> nf2 = dynamic_pointer_cast<nonlinear_function>(f2);

	// Check
	if(!nf1 || !nf2) {
		std::cerr << "<<ERROR>> One of the input function of the product is NULL" << std::endl;
		return ptr<function>(NULL);
	}

	product_function* pf = new product_function(nf1, nf2);
	return ptr<function>(pf);
}

/* Setting/Get the parameters of a function object
 *
 * TODO: Add the rational function interface
 */
void set_function_params(ptr<function>& f, const vec& v) {

	// Try to set the parameter as a nonlinear function
	ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(f);
	if(nf) {
		if(nf->nbParameters() == v.size()) {
			nf->setParameters(v);
		} else {
			std::cerr << "<<ERROR>> Vector of params and function have different sizes." << std::endl;
		}
		return;
	}

	// ptr<rational_function> rf = dynamic_pointer_cast<rational_function>(f);
	// if(rf) {
	// 	int np, nq;
	// 	rf->size(np, nq);

	// 	if(np+nq == v.size()) {
	// 		rf->setParameters(v);
	// 	} else {
	// 		std::cerr << "<<ERROR>> Vector of parameters has different size that the functions number of parameters" << std::endl;
	// 	}
	// 	return;
	// }
}

python_vec get_function_params(ptr<function>& f) {
		// Try to set the parameter as a nonlinear function
	ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(f);
	if(nf) {
		return nf->parameters();
	}

	std::cerr << "<<ERROR>> Parameters cannot be retrieved" << std::endl;
	vec res(1);
	return res;
}

/* Save a function object to a file, without any argument option. This will save the function
 * object in ALTA's format.
 */
void save_function_without_args(const ptr<function>& f, const std::string& filename) {
   arguments args;
   f->save(filename, args);
}

/* Fitter interface to allow to launch fit with and without providing an
 * arguments object.
 */
bool fit_data_without_args(const ptr<fitter>& _fitter, const ptr<data>& _data, ptr<function>& _func) {
   arguments args;
   _fitter->set_parameters(args);
   return _fitter->fit_data(_data, _func, args);
}

bool fit_data_with_args(ptr<fitter>& _fitter, const ptr<data>& _data, ptr<function>& _func, const arguments& args) {
   _fitter->set_parameters(args);
   return _fitter->fit_data(_data, _func, args);
}


/* Softs functions. Those function recopy the softs main function, without
 * the command line arguments.
 * TODO: Add the command line arguments in the parameters
 */
void data2data(const data* d_in, data* d_out)
{
   if(dynamic_cast<vertical_segment*>(d_out)!=NULL)
   {
      d_out->setParametrization(d_in->input_parametrization());
      d_out->setDimX(d_in->dimX());
      d_out->setDimY(d_in->dimY());

      // Init the min and max
      vec _min(d_out->dimX()), _max(d_out->dimX());
      for(unsigned int k=0; k<d_out->dimX(); ++k)
      {
         _min[k] =  std::numeric_limits<double>::max() ;
         _max[k] = -std::numeric_limits<double>::max() ;
      }

      vec temp(d_out->dimX() + d_out->dimY());
      for(unsigned int i=0; i<d_in->size(); ++i)
      {
         // Copy the input vector
         vec x = d_in->get(i);
         params::convert(&x[0], d_in->parametrization(), d_out->parametrization(), &temp[0]);
         params::convert(&x[d_in->dimX()], d_in->output_parametrization(), d_in->dimY(), d_out->output_parametrization(), d_out->dimY(), &temp[d_out->dimX()]);
         d_out->set(temp);

         // Update min and max
         for(unsigned int k=0; k<d_out->dimX(); ++k)
         {
            _min[k] = std::min(_min[k], temp[k]) ;
            _max[k] = std::max(_max[k], temp[k]) ;
         }
      }

      // Set the min and max
      d_out->setMin(_min);
      d_out->setMax(_max);
   }
   else
   {
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
}

/* This function provides a similar behaviour that the brdf2data function.
 * An input function object is evaluated on the input position of a data
 * object. To do so, the data object must contains some positions. It can
 * be so when the data object is a laoded data sample or when the data type
 * has predefined sample sets.
 */
void brdf2data(const ptr<function>& f, ptr<data>& d) {
	if(d->size() == 0) {
		std::cerr << "<<ERROR>> Please provide a data object with a sample structure or load a data file with defined positions." << std::endl;
		return;
	}

	vec temp(f->dimX());
	for(int i=0; i<d->size(); ++i) {
		// Copy the input vector
		vec x = d->get(i);

		// Convert the data to the function's input space.
	    if(f->input_parametrization() == params::UNKNOWN_INPUT) {
	    	memcpy(&temp[0], &x[0], f->dimX()*sizeof(double));
	    } else {
			params::convert(&x[0], d->parametrization(), f->parametrization(), &temp[0]);
		}
		vec y = f->value(temp);

		for(int j=0; j<d->dimY(); ++j) {
			x[d->dimX() + j] = y[j];
		}

		d->set(i, y);
	}	
}


/*! \inpage python 
 *  Exporting the ALTA module 
 */
BOOST_PYTHON_MODULE(alta)
{
	// Custor converters
	//
	list_to_vec_converter();
	vec_to_list_converter();


	// Argument class
	//
	bp::class_<arguments>("_arguments");
	bp::class_<python_arguments, bp::bases<arguments>>("arguments")
		.def(bp::init<>())
		.def(bp::init<bp::dict>())
		.def("update", &arguments::update);


	// Vector class
	//
	// TODO: There is a conversion issue right now that prevents us from using vectors
	// within Python. This needs to be investiguated.
	bp::class_<vec>("_vec")
		.def("__len__", &vec::size)
		.def("__getitem__", &vec_get_item)
		.def("__setitem__", &vec_set_item);
	bp::class_<python_vec>("vec")
		.def(bp::init<vec>())
		.def(bp::init<bp::list>())
		.def(bp::self_ns::str(bp::self_ns::self));
	bp::implicitly_convertible<vec, python_vec>();
	bp::implicitly_convertible<python_vec, vec>();


	// Function interface
	//
	bp::class_<function, ptr<function>, boost::noncopyable>("function", bp::no_init)
		.def("__add__", &add_function)
		.def("__mul__", &mult_function)
		.def("__rmul__", &mult_function)
		.def("value", &function::value)
		.def("load", &load_from_file)
		.def("load", &load_from_file_with_args)
		.def("save",  &function::save)
      .def("save", &save_function_without_args)
		.def("set", &set_function_params)
		.def("get", &get_function_params);
	bp::def("get_function", get_function);
	bp::def("get_function", get_function_from_args);
	bp::def("load_function", load_function);
	bp::def("load_function", load_function_with_args);



	// Data interface
	//
	bp::class_<data, ptr<data>, boost::noncopyable>("data", bp::no_init)
		.def("size", &data::size)
		.def("get",  &data::get)
		//.def("set",  &data::set)
		.def("load", static_cast< void(data::*)(const std::string&)>(&data::load))
		.def("load", static_cast< void(data::*)(const std::string&, const arguments&)>(&data::load))
		.def("save", &data::save);
	bp::def("get_data",  get_data);
	bp::def("get_data",  get_data_with_args);
	bp::def("load_data", load_data);


	// Fitter interface
	//
	bp::class_<fitter, ptr<fitter>, boost::noncopyable>("fitter", bp::no_init)
		.def("fit_data", &fit_data_with_args)
      .def("fit_data", &fit_data_without_args);
	bp::def("get_fitter", plugins_manager::get_fitter);

	// Softs
	//
	bp::def("data2data", data2data);
	bp::def("brdf2data", brdf2data);	
}
