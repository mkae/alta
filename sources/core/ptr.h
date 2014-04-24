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

/*! \brief Define a shared pointer class.
 */
template<class T> class ptr
{
	public:
		//! Empty constructor
		ptr(T* ptr) : _ptr(ptr)
	   {
			_counter = new ptr_counter();
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
				delete _ptr;
				delete _counter;
			}
		}

		//! Evaluation operator. This operator should be inlined for 
		//! performance reasons.
		inline T* operator-> ()
		{
			return _ptr;
		}

	private:
		T* _ptr;
		ptr_counter* _counter;
};
