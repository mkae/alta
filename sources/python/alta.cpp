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


// here comes the magic
template <typename T> 
T* get_pointer(ptr<T> const& p) {
  //notice the const_cast<> at this point
  //for some unknown reason, bp likes to have it like that
  return const_cast<T*>(p.get());
}

// some boost.python plumbing is required as you already know
namespace boost { namespace python {

    template <typename T>
    struct pointee< ::ptr<T> > {
        typedef T type;
    };

}}

struct my_vec : public vec {
    my_vec() : vec() {}

    my_vec(const bp::list& l) : vec(bp::len(l)) {
        for(size_t i=0; i<bp::len(l); ++i) {
            (*this)[i] = bp::extract<double>(l[i]);
        }
    }
};

ptr<data> load_data(const std::string& plugin_name, const std::string& filename) {
	ptr<data> d = plugins_manager::get_data(plugin_name);
	d->load(filename);
	return d;
}

// Creating functions for the plugins_manager calls
function* get_function(const std::string& filename) {
    return plugins_manager::get_function(filename);
}
function* get_function_from_args(const arguments& args) {
    return plugins_manager::get_function(args);
}

/* Exporting the ALTA module */
BOOST_PYTHON_MODULE(alta)
{
	// Argument class
	//
	// TODO: {Laurent: we should be able to set arguments using python's maps}
	bp::class_<arguments>("arguments")
		.def(bp::init<>())
		.def("update", &arguments::update);

	bp::class_<my_vec>("vec")
		.def(bp::init<bp::list>())
		/*        .def(" __setitem__", &my_vec::operator[])*/;

	// Function interface
	//
	bp::class_<function, ptr<function>, boost::noncopyable>("function", bp::no_init)
		.def("value", &function::value);
	bp::def("get_function", get_function, bp::return_value_policy<bp::manage_new_object>());
	bp::def("get_function", get_function_from_args, bp::return_value_policy<bp::manage_new_object>());

	// Data interface
	//
	bp::class_<data, ptr<data>, boost::noncopyable>("data", bp::no_init)
		.def("load", static_cast< void(data::*)(const std::string&)>(&data::load))
		.def("size", &data::size);
	bp::def("get_data", plugins_manager::get_data);
	bp::def("load_data", load_data);

	// Fitter interface
	//
	// TODO: {Laurent: to make the call to fit_data possible, we need to
	//        implement a Python converter to ptr<T>. This seems to be
	//        tricky (see boost/python/shared_ptr_to_python.hpp). The
	//        other solution is to use boost's shared_ptr instead of ptr
	//        when compiling for the Python interface.}
	bp::class_<fitter, ptr<fitter>, boost::noncopyable>("fitter", bp::no_init)
		.def("fit_data", &fitter::fit_data);
	bp::def("get_fitter", plugins_manager::get_fitter);
}
