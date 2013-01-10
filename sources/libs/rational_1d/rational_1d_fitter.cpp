#include "rational_1d_fitter.h"

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>
#include <iostream>
#include <fstream>

rational_1d::rational_1d() 
{
}

rational_1d::rational_1d(const std::vector<float>& a,
                         const std::vector<float>& b) :
	a(a), b(b)
{
}
rational_1d::~rational_1d()
{
}

// Overload the function operator
float rational_1d::operator()(float x) const 
{
	float p = 0.0f ;
	float q = 0.0f ;

	for(int i=a.size()-1; i>=0; --i)
	{
		p = x*p + a[i] ;
	}

	for(int i=b.size()-1; i>=0; --i)
	{
		q = x*q + b[i] ;
	}

	return p/q ;
}
		
// Get the p_i and q_j function
float rational_1d::p(float x, int i) const
{
	return pow(x, i) ;
}
float rational_1d::q(float x, int j) const 
{
	return pow(x, j) ;
}

// IO function to text files
void rational_1d::load(const std::string& filename)
{
}
void rational_1d::save() const
{
}

void rational_1d_data::load(const std::string& filename) 
{
	std::ifstream file(filename) ;

	float x, y, dy ;
	while(file.good())
	{
		file >> x >> y >> dy ;
		std::vector<float> v ;
		v.push_back(x) ;
		v.push_back(y-dy) ;
		v.push_back(y+dy) ;
		data.push_back(v) ;
	}
}

bool rational_1d_data::get(int i, float& x, float& yl, float& yu) const
{
	if(i >= data.size())
	{
		return false ;
	}

	x  = data[i][0] ;
	yl = data[i][1] ;
	yu = data[i][2] ;

	return true ;
}

const std::vector<float>& rational_1d_data::operator[](int i) const
{
	return data[i] ;
}

int rational_1d_data::size() const
{
	return data.size() ;
}

typedef CGAL::MP_Float ET ;
typedef CGAL::Quadratic_program<float> Program ;
typedef CGAL::Quadratic_program_solution<ET> Solution ;
		
// Fitting a data
rational_1d rational_1d_fitter::fit_data(const rational_1d_data& data)
{
	// The resulting ration function
	rational_1d r ;

	// by default, we have a nonnegative QP with Ax <= b
	Program qp (CGAL::LARGER, false, 0, false, 0) ; 

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	int np = 10, nq = 10;
	for(int i=0; i<np+nq; ++i)
	{
		qp.set_d(i, i, 1.0f) ; 
	}

	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	for(int i=0; i<data.size(); ++i)	
	{		
		// Norm of the row vector
		float a0_norm = 0.0f ;
		float a1_norm = 0.0f ;

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<np+nq; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const float pi = r.p(data[i][0], j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				qp.set_a(j, 2*i+0,  pi) ;
				qp.set_a(j, 2*i+1, -pi) ;
			}
			// Filling the q paret
			else
			{
				const float qi = r.q(data[i][0], j-np) ;
				a0_norm += qi*qi * (data[i][1]*data[i][1]) ;
				qp.set_a(j, 2*i+0, -data[i][1] * qi) ;
				
				a1_norm += qi*qi * (data[i][2]*data[i][2]) ;
				qp.set_a(j, 2*i+1,  data[i][2] * qi) ;
			}
		}
	
		// Set the c vector, in Oliver's work it is
		// using the max SVD value.
		qp.set_b(2*i+0, a0_norm) ;
		qp.set_b(2*i+1, a1_norm) ;
	}

#ifdef DEBUG
	// Export some informations on the problem to solve
	std::cout << qp.get_n() << " variables" << std::endl ;
	std::cout << qp.get_m() << " constraints" << std::endl ;
#endif

	if(qp.is_linear() && !qp.is_nonnegative())
	{
		std::cerr << "Error here, should not be linear or negative!" << std::endl ;
		throw ;
	}

#ifdef DEBUG
	// Iterate over the rows
	std::cout << "A = [" ;
	for(int j=0; j<qp.get_m(); ++j)
	{
		if(j == 0) std::cout << "   " ;

		// Iterate over the columns
		for(int i=0; i<qp.get_n(); ++i)
		{
			if(i > 0) std::cout << ",\t" ;
			std::cout << *(*(qp.get_a()+i) +j) ;
		}
		std::cout << ";" ;
		if(j < qp.get_n()-1) std::cout << std::endl ;
	}
	std::cout << "]" << std::endl ;
#endif

	// solve the program, using ET as the exact type
	Solution s = CGAL::solve_quadratic_program(qp, ET()) ;
	assert (s.solves_quadratic_program(qp)) ;

	// Recopy the vector data
	std::vector<float> p, q;
	for(int i=0; i<np+nq; ++i)
	{
		const float v = (float)CGAL::to_double(*(s.variable_numerators_begin()+i)) ;

		if(i < np)
		{
			p.push_back(v) ;
		}
		else
		{
			q.push_back(v) ;
		}
	}

	r = rational_1d(p, q);
	return r;
}
