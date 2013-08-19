template <class RF1D>
rational_function_t<RF1D>::rational_function_t() : np(0), nq(0)
{
}

template <class RF1D>
rational_function_t<RF1D>::rational_function_t(int np, int nq) : np(np), nq(nq)
{
}

//! \todo clean memory here
template <class RF1D>
rational_function_t<RF1D>::~rational_function_t()
{
}

template <class RF1D>
void rational_function_t<RF1D>::update(int i,
                                       RF1D* r)
{
    rs[i] = r;
}

template <class RF1D>
RF1D* rational_function_t<RF1D>::get(int i)
{
    // Check for consistency in the index of color channel
    if(i < _nY)
    {
        if(rs[i] == NULL)
        {
            rs[i] = new RF1D(np, nq);
            rs[i]->setDimX(dimX());
            rs[i]->setDimY(dimY());
            rs[i]->setMin(getMin()) ;
            rs[i]->setMax(getMax()) ;
        }
        return rs[i];
    }
    else
    {
        std::cout << "<<ERROR>> tried to access out of bound 1D RF" << std::endl;
        return NULL;
    }
}

// Overload the function operator
template <class RF1D>
vec rational_function_t<RF1D>::value(const vec& x) const
{
    vec res(_nY) ;

    for(int k=0; k<_nY; ++k)
    {
        res[k] = rs[k]->value(x)[0] ;
    }
    return res ;
}

// IO function to text files
template <class RF1D>
void rational_function_t<RF1D>::load(const std::string& filename)
{
    std::ifstream file(filename.c_str()) ;
    if(!file.is_open())
    {
        std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
        throw ;
    }

    int nX, nY ;
    vec xmin, xmax ;
    int i = 0, j = 0;
    while(file.good())
    {
        std::string line ;
        std::getline(file, line) ;
        std::stringstream linestream(line) ;

        // Discard incorrect lines
        if(linestream.peek() == '#')
        {
            linestream.ignore(1) ;

            std::string comment ;
            linestream >> comment ;

            if(comment == std::string("DIM"))
            {
                linestream >> nX >> nY ;
                setDimX(nX) ;
                setDimY(nY) ;

                xmin.resize(nX) ;
                xmax.resize(nX) ;
                for(int k=0; k<nX; ++k)
                    xmax[k] = 1.0;

                setMin(xmin) ;
                setMax(xmax) ;
            }
            else if(comment == std::string("NP"))
            {
                linestream >> np ;
                //a.resize(np*nY);
            }
            else if(comment == std::string("NQ"))
            {
                linestream >> nq ;
                //b.resize(nq*nY);
            }
            else if(comment == std::string("INPUT_PARAM"))
            {
                std::string param;
                linestream >> param ;


                setParametrization(params::parse_input(param));
            }
            continue ;
        }
        else if(line.empty())
        {
            continue ;
        }
        else if(j < nY)
        {
            int index ; double val ;

            // Accessing the index
            for(int k=0; k<nX; ++k) {
                linestream >> index ;
            }

            // Accessing the value
            linestream >> val ;
            /*! \todo
            if(i < np)
            {
                a[i + np*j] = val ;
            }
            else
            {
                b[i-np + nq*j] = val ;
            }

            if(i < np+nq) {
                ++i ;
            } else {
                i = 0 ;
                ++j ;
            }
            */
        }
    }
/*
    for(int i=0; i<a.size(); ++i) {
        std::cout << a[i] << "\t" ;
    }
    for(int i=0; i<b.size(); ++i) {
        std::cout << b[i] << "\t" ;
    }
*/
}

//! \todo it should handle parametrization
template <class RF1D>
void rational_function_t<RF1D>::save_matlab(const std::string& filename, const arguments& args) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);

    file << "function y = brdf(x)" << std::endl;
    file << std::endl;
    file << "\ts = [";
    for(int i=0; i<dimX(); ++i)
    {
        file << 1.0 / (_max[i]-_min[i]);
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "];" << std::endl;
    file << "\tc = [";
    for(int i=0; i<dimX(); ++i)
    {
        file << _min[i];
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "];" << std::endl;
    file << std::endl ;

    // Export each color channel independantly
    for(int j=0; j<dimY(); ++j)
    {
        // Export the numerator of the jth color channel
        file << "\tp(" << j+1 << ",:) = ";
        for(unsigned int i=0; i<np; ++i)
        {
              /*! \todo
            if(i > 0 && a[np*j + i] >= 0.0)
                file << " + ";
            else if(a[np*j + i] < 0.0)
                file << " " ;

            file << a[np*j + i];

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
               file << ".*legendrepoly(" << degree[k] << ", 2.0*((x(" << k+1 << ",:)"
                         << "-c(" << k+1 << "))*s(" << k+1 << ") - 0.5))" ;
            }
                */
        }
        file << ";" << std::endl;

        // Export the denominator of the jth color channel
        file << "\tq(" << j+1 << ",:) = ";
        for(unsigned int i=0; i<nq; ++i)
        {
              /*! \todo
            if(i > 0 && b[np*j + i] >= 0.0)
                file << " + ";
            else if(b[np*j + i] < 0.0)
                file << " " ;

            file << b[np*j + i] ;

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
               file << ".*legendrepoly(" << degree[k] << ", 2.0*((x(" << k+1 << ",:)"
                         << "-c(" << k+1 << "))*s(" << k+1 << ") - 0.5))" ;
            }
                */
        }
        file << ";" << std::endl;

        file << "\ty(" << j+1 << ",:) = p./q;" << std::endl;
        if(j < dimY()-1)
        {
            file << std::endl;
        }
    }


    file << "endfunction" << std::endl;

    file.close() ;
}

//! \todo it should handle parametrization
template <class RF1D>
void rational_function_t<RF1D>::save_cpp(const std::string& filename, const arguments& args) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);

    file << "double s[" << dimX() << "] = {";
    for(int i=0; i<dimX(); ++i)
    {
        file << 1.0 / (_max[i]-_min[i]);
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "};" << std::endl;
    file << "double c[" << dimX() << "] = {";
    for(int i=0; i<dimX(); ++i)
    {
        file << _min[i];
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "};" << std::endl;
    file << std::endl ;

    file << "// The Legendre polynomial of order i evaluated in x" << std::endl;
    file << "double l(double x, int i)" << std::endl;
    file << "{" << std::endl;
    file << "    if(i == 0)" << std::endl;
    file << "    {" << std::endl;
    file << "        return 1;" << std::endl;
    file << "    }" << std::endl;
    file << "    else if(i == 1)" << std::endl;
    file << "    {" << std::endl;
    file << "        return x;" << std::endl;
    file << "    }" << std::endl;
    file << "    else" << std::endl;
    file << "    {" << std::endl;
    file << "        return ((2*i-1)*x*l(x, i-1) - (i-1)*l(x, i-2)) / (double)i ;" << std::endl;
    file << "    }" << std::endl;
    file << "}" << std::endl;
    file << std::endl;

    file << "void brdf(double* x, double* y)" << std::endl;
    file << "{" << std::endl;
    file << "\tdouble p, q;" << std::endl;

    // Export each color channel independantly
    for(int j=0; j<dimY(); ++j)
    {
        // Export the numerator of the jth color channel
        file << "\tp = ";
        for(unsigned int i=0; i<np; ++i)
        {
              /*! \todo
            if(i > 0 && a[np*j + i] >= 0.0)
            {
                file << " + ";
            }
            file << a[np*j + i];

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
                file << "*l(2.0*((x[" << k << "]-c[" << k << "])*s[" << k << "] - 0.5), " << degree[k] << ")" ;
            }
                */
        }
        file << ";" << std::endl;

        // Export the denominator of the jth color channel
        file << "\tq = ";
        for(unsigned int i=0; i<nq; ++i)
        {
              /*! \todo
            if(i > 0 && b[np*j + i] >= 0.0)
                file << " + ";
            else if(b[np*j + i] < 0.0)
                file << " " ;

            file << b[np*j + i] ;

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
                file << "*l(2.0*((x[" << k << "]-c[" << k << "])*s[" << k << "] - 0.5), " << degree[k] << ")" ;
            }
                */
        }
        file << ";" << std::endl;

        file << "\ty[" << j << "] = p/q;" << std::endl;
        if(j < dimY()-1)
        {
            file << std::endl;
        }
    }


    file << "}" << std::endl;

    file.close() ;
}

template <class RF1D>
void rational_function_t<RF1D>::save_gnuplot(const std::string& filename, const data* d, const arguments& args) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);
    for(int i=0; i<d->size(); ++i)
    {
        vec v = d->get(i) ;
//		vec y1 ; y1.assign(d->dimY(), 0.0) ;
//		for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

        vec y2 = value(v) ;
        for(int u=0; u<d->dimX(); ++u)
            file << v[u] << "\t" ;

        for(int u=0; u<d->dimY(); ++u)
            file << y2[u] << "\t" ;

        file << std::endl ;
    }
    file.close();
}

template <class RF1D>
void rational_function_t<RF1D>::save(const std::string& filename) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);
    file << "#DIM " << _nX << " " << _nY << std::endl ;
    file << "#NP " << np << std::endl ;
    file << "#NQ " << nq << std::endl ;
    file << "#BASIS LEGENDRE" << std::endl ;
    file << "#INPUT_PARAM " << params::get_name(this->parametrization()) << std::endl;

    for(int k=0; k<_nY; ++k)
    {
        /*
        for(unsigned int i=0; i<np; ++i)
        {
            std::vector<int> index = index2degree(i) ;
            for(unsigned int j=0; j<index.size(); ++j)
            {
                file << index[j] << "\t" ;
            }
            file << a[i+np*k] << std::endl ;
        }
        */

        /*
        for(unsigned int i=0; i<nq; ++i)
        {
            std::vector<int> index = index2degree(i) ;
            for(unsigned int j=0; j<index.size(); ++j)
            {
                file << index[j] << "\t" ;
            }
            file << b[i+nq*k] << std::endl ;
        }
        */
    }

}

template <class RF1D>
std::ostream& operator<< (std::ostream& out, rational_function_t<RF1D>& r)
{
    for(int i=0; i<r.dimY(); ++i)
    {
        RF1D* rf = r.get(i);
        out << "dimension " << i << ": ";
        if(rf != NULL)
        {
            out << *rf << std::endl;
        }
        else
        {
            out << "[NULL]" << std::endl;
        }
    }

    return out ;
}
