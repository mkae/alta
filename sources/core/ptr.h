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
#if __cplusplus >= 201103L

#include <memory>
template<class T> using ptr = std::shared_ptr<T>;

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
			if(_counter->value() == 0)
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
