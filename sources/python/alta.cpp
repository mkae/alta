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
template <typename T> T* get_pointer(ptr<T> const& p) {
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

// Creating a non purely virtual class for the function
struct Function : function, bp::wrapper<function> { 
    vec value(const vec& x) const {
        return this->get_override("value")(x);
    }

    bool load(std::istream& in) {
        return this->get_override("load")(in);
    }
};

// Creating a non purely virtual class for the data
struct Data : data, bp::wrapper<data> { 
		// Load data from a file
		void load(const std::string& filename) {
        this->get_override("load")(filename);
		}
		void load(const std::string& filename, const arguments& args) {
        this->get_override("load")(filename, args);
		}

		vec get(int i) const {
        return this->get_override("get")(i);
		}

		vec operator[](int i) const {
        return this->get_override("operator[]")(i);
		}
		
		vec value(vec in) const {
        return this->get_override("value")(in);
		}

		void set(vec x) {
        this->get_override("set")(x);
		}

		int size() const {
        return this->get_override("size")();
		}
};

ptr<data> load_data(const std::string& plugin_name, const std::string& filename) {
	ptr<data> d = plugins_manager::get_data(plugin_name);
	d->load(filename);
	return d;
}

// Creating a non purely virtual class for the function
struct Fitter : fitter, bp::wrapper<fitter> { 
	bool fit_data(const ptr<data>& d, ptr<function>& f, const arguments& args) {
		return this->get_override("fit_data")(f, args);
	}

	void set_parameters(const arguments& args) {
		this->get_override("set_parameters")(args);
	}
};

// Creating functions for the plugins_manager calls
function* get_function(const std::string& filename) {
    return plugins_manager::get_function(filename);
}

// Exporting the ALTA module
BOOST_PYTHON_MODULE(alta)
{
    bp::class_<arguments>("arguments")
		 .def(bp::init<>());

	 bp::class_<my_vec>("vec")
        .def(bp::init<bp::list>())
/*        .def(" __setitem__", &my_vec::operator[])*/;

    bp::class_<Function, ptr<Function>, boost::noncopyable>("function")
    	.def("value", bp::pure_virtual(&function::value));

	 // Function interface
    bp::def("get_function", get_function, bp::return_value_policy<bp::manage_new_object>());

	 // Data interface
    bp::class_<Data, ptr<Data>, boost::noncopyable>("data")
//		 .def("load", bp::pure_virtual(&data::load));
		 .def("size", bp::pure_virtual(&data::size));
	 bp::def("get_data", plugins_manager::get_data);
	 bp::def("load_data", load_data);
	 bp::register_ptr_to_python< ptr<data> >();

	 // Fitter interface
    bp::class_<Fitter, ptr<Fitter>, boost::noncopyable>("fitter")
		 .def("fit_data", bp::pure_virtual(&fitter::fit_data));
	 bp::def("get_fitter", plugins_manager::get_fitter);
	 bp::register_ptr_to_python< ptr<fitter> >();
}
