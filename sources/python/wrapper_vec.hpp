/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Université de Montréal

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// Boost includes
#include <boost/python.hpp>

// ALTA include
#include <core/common.h>

// STL include
#include <iostream>

#define bp boost::python

/* Create special converter from and to list/vec
 *
 * This allow to create a 'vec' object using a python 'list' with the following
 * interface:
 *
 *    v = alta.vec([1,1,1])
 *
 */
struct converter_vec_list {

	// Constructor add the converter to the registry
	static void register_converter() {
		bp::converter::registry::push_back(
				&convertible,
				&construct,
				bp::type_id<vec>());
	}

	// Is the Python object convertible to a vec ?
	static void* convertible(PyObject* obj_ptr) {
		if (!PyList_Check(obj_ptr)) return nullptr;
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

/* Iterator over the 'vec' class to enable construction of lists and arrays in
 * python using the __iter__ interface. To build a 'list' or a numpy 'array'
 * from an object (like 'vec'), Python requires another class that enable to
 * iterate over the structure. This class is created using 'vec::__iter__'
 * method. Then, elements of 'vec' are extracted using the 'next' method of the
 * iterator. To enable use with numpy.array, the iterator needs to be nested
 * (have a method '__iter__'). This is how numpy handles multidmensional arrays.
 *
 * With this class, it is possible to use the following code:
 *
 *    v = alta.vec([1,1,1])
 *    a = np.array(v)
 */
struct iterator_vec {
   const double* cur;
   const double* end;

   iterator_vec(const vec& x) : cur(x.data()), end(x.data()+x.size()) { }

   // Iterator over the 'vec' structure.
   //
   static double next(iterator_vec& iter) {
      if(iter.cur == iter.end) {
         PyErr_SetString(PyExc_StopIteration, "No more data.");
         bp::throw_error_already_set();
      }

      double r = *(iter.cur);
      iter.cur++;
      return r;
   }

   // Create an iterator from the 'vec' object.
   //
   static iterator_vec iter(const vec& x) {
      return iterator_vec(x);
   }

   // Create a nested iterator.
   //
   static iterator_vec self(const iterator_vec& x) {
      return x;
   }
};

/* Specific function call to acces a vector's element
 */
double vec_get_item(const vec& x, int i) {
	return x[i];
}

/* Specific function call to set a vector's element
 */
void vec_set_item(vec& x, int i, double a) {
	x[i] = a;
}

/* Operators on vec
 */
inline vec vec_add(const vec& a, const vec& b) {
   return a + b;
}
inline vec vec_sub(const vec& a, const vec& b) {
   return a - b;
}

/* Specific convert a vec to a string
 */
std::string vec_str(const vec& x) {
   std::stringstream stream;
   stream << "[";
   for(int i=0; i<x.size(); ++i) {
      if(i > 0) { stream << ", "; }
      stream << x[i];
   }
   stream << "]";
   return stream.str();
}

/* Register the vector class to boost python
 */
void register_wrapper_vec() {

   converter_vec_list::register_converter();

	bp::class_<iterator_vec>("iterator_vec", bp::no_init)
      .def("next", &iterator_vec::next)
      .def("__iter__", &iterator_vec::self);

	bp::class_<vec>("vec")
		.def(bp::init<vec>())
		.def(bp::self_ns::str(bp::self_ns::self))
		.def("__add__", &vec_add)
		.def("__sub__", &vec_sub)
		.def("__len__", &vec::size)
		.def("__getitem__", &vec_get_item)
		.def("__setitem__", &vec_set_item)
		.def("__str__", &vec_str)
      .def("__iter__", &iterator_vec::iter);
}
