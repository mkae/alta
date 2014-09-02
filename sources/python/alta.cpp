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
struct Function : function, bp::wrapper<function>
{ 
    vec value(const vec& x) const
    {
        return this->get_override("value")(x);
    }

    bool load(std::istream& in) {
        return this->get_override("load")(in);
    }
};

// Creating functions for the plugins_manager calls
function* get_function(const std::string& filename) {
    return plugins_manager::get_function(filename);
}

// Exporting the ALTA module
BOOST_PYTHON_MODULE(alta)
{
    bp::class_<my_vec>("vec")
        .def(bp::init<bp::list>())
/*        .def(" __setitem__", &my_vec::operator[])*/;

    bp::class_<Function, ptr<Function>, boost::noncopyable>("function")
    	.def("value", bp::pure_virtual(&function::value));

    bp::def("get_function", get_function, bp::return_value_policy<bp::manage_new_object>());

}
