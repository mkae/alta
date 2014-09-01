// Boost includes
#include <boost/python.hpp>

// ALTA include
#include <core/function.h>

using namespace boost::python;


// Creating a non purely virtual class for the function
struct Function : function, wrapper<function>
{ 
    vec value(const vec& x) const
    {
        return this->get_override("value")(x);
    }

    bool load(std::istream& in) {
        return this->get_override("load")(in);
    }
};

// Exporting the ALTA module
BOOST_PYTHON_MODULE(alta)
{

    class_<Function, boost::noncopyable>("function")
    	.def("value", pure_virtual(&function::value))
    ;


}
