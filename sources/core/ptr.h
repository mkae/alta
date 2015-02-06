/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

/*! \class ptr
 *  
 *  Define a shared pointer class for automatic memory cleaning
 *  when dealing with objects created by plugins. All plugins should
 *  provide a ptr<object> instead of a object* so that the user won't
 *  have to worry about allocation/deallocation.
 *
 *  ALTA provide a weak implementation of the ptr class. If a C++11
 *  compatible is used, it will try to use the STL version of shared
 *  pointer.
 */

/* Checking for the presence of C++11 features like smart pointers. There
 * is no clean way to do it for all compilers. This method is supposed to
 * work with CLANG and GCC.
 *
 * See htpp://stackoverflow.com/questions/11886288/ 
 */
#if defined(USE_BOOST)

#include <boost/shared_ptr.hpp>
#define ptr boost::shared_ptr
#define dynamic_pointer_cast boost::dynamic_pointer_cast

#elif defined(TODO)
 //__cplusplus >= 201103L

#include <memory>
template<class T> using ptr = std::shared_ptr<T>;

template<class T, class U> 
inline ptr<U> dynamic_pointer_cast(const ptr<T>& ptr_t) {
	return std::dynamic_pointer_cast<U>(ptr_t);
}

#else

/*  Define a counter class. This class should not be used by any other part of
 *  ALTA code.
 */
struct ptr_counter
{
	ptr_counter() : _count(1) { }

	inline void increment()
	{
		++_count;
	}

	inline void decrement()
	{
		--_count;
	}

	inline unsigned int value() const
	{
		return _count;
	}

	unsigned int _count;
};

template<class T> class ptr;

template<class T, class U> 
ptr<U> dynamic_pointer_cast(const ptr<T>& ptr_t);

template<class T> class ptr
{
	public:
		//! Empty constructor
		ptr(T* ptr) : _ptr(ptr)
	   {
			_counter = new ptr_counter();
		}

		//! Counter copy constructor
		ptr(T* ptr, ptr_counter* counter) : _ptr(ptr), _counter(counter)
		{
			_counter->increment();	
		}

		//! Copy constructor
		ptr(const ptr<T>& p) : _ptr(p._ptr), _counter(p._counter)
		{
			_counter->increment();
		}

		//! Destructor, should free memory when the counter goes
		//! to zero.
		~ptr()
		{
			_counter->decrement();
			if(_counter->value() < 1 || _ptr == NULL)
			{
				if(_ptr != NULL) { delete _ptr; }
				delete _counter;
			}
		}

		//! Evaluation operator. This operator should be inlined for 
		//! performance reasons.
		inline T* operator-> () const
		{
			return _ptr;
		}

		//! Assignment operator. If a valid pointer is already present, its
		//! counter is decremented and the pointer and counter might be
		//! deleted in case the counter reach zero. Then, the elements of
		//! a, both counter and pointer, are copied to this and the counter
		//! is incremented.
		ptr<T>& operator=(const ptr<T>& a)
		{
			//RP:  Check and avoid assignment to itself 
		    if((void*) this == (void *) &a)
		    {
		        return *this;
		    }

			_counter->decrement();
			if(_counter->value() < 1)
			{
				if(_ptr != NULL) { delete _ptr; }
				delete _counter;
			}

			_counter = a._counter;
			_ptr     = a._ptr;
			_counter->increment();

			return *this;
		}

		//! Raw acces to the pointer. It is sometimes needed to 
		inline T* get() const
		{
			return _ptr;
		}

		//! Is the underlying pointer not NULL.
		inline operator bool() const
		{
			return _ptr != NULL;
		}

		template<class U> 
		friend ptr<U> dynamic_pointer_cast(const ptr<T>& ptr_t)
		{
			U* u = dynamic_cast<U*>(ptr_t._ptr);
			if(u == NULL)
			{
				return ptr<U>(NULL);
			}
			else
			{
				ptr<U> ptr_u = ptr<U>(u, ptr_t._counter);
				return ptr_u;
			}
		}

	private:

		T* _ptr;
		ptr_counter* _counter;

};




#endif
