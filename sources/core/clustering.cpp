#include "clustering.h"

#include <iostream>
#include <limits>
#include <string>
#include <fstream>

clustering::clustering(const data* d, const arguments& args)
{
    // List of the indices of the clustering
    std::vector<int> indices;

    vec findices = args.get_vec("cluster-dim", d->dimX(), -1.0f);
    for(int i=0; i<findices.size() && findices[i] != -1.0f; ++i)
    {
        int id = (int)findices[i];
        indices.push_back(id);

        assert(id >= 0 && id < d->dimX());
    }
    std::cout << "<<DEBUG>> the resulting data will have " << indices.size() << " dimensions" << std::endl;

    // We cannot have null cluster size of full rank.
    assert(indices.size() > 0 && indices.size() < d->dimX());

    // Fit the slice into the domain of definition of the data
    vec compressed_vec  = args.get_vec("cluster-slice", d->dimX(), std::numeric_limits<float>::max());
    for(int i=0; i<d->dimX(); ++i)
    {
        compressed_vec[i] = std::min<double>(compressed_vec[i], d->max()[i]);
        compressed_vec[i] = std::max<double>(compressed_vec[i], d->min()[i]);
    }

    // Set the limit of the new domain and dimensions
    _nX = indices.size();
    _nY = d->dimY();
    _min = vec(_nX);
    _max = vec(_nX);
    for(int i=0; i<_nX; ++i)
    {
        _min[i] = d->min()[indices[i]];
        _max[i] = d->max()[indices[i]];
    }

    for(int i=0; i<d->size(); ++i)
    {
        vec p = d->get(i);
        bool reject = false;
        for(int i=0; i<d->dimX(); ++i)
        {
            if(is_in<int>(indices, i) == -1 && p[i] != compressed_vec[i])
            {
                reject = true;
            }
        }
        if(reject)
        {
            continue;
        }

        vec e(dimX()+dimY());
        for(int i=0; i<_nX; ++i)
        {
            e[i] = p[indices[i]];
        }
        for(int i=0; i<_nY; ++i)
        {
            e[_nX + i] = p[d->dimX() + i];
        }
        _data.push_back(e);
    }
    std::cout << "<<INFO>> clustering left " << _data.size() << " elements" << std::endl;

    save("cluster.txt");
}

void clustering::save(const std::string& filename)
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);
    file << "#DIM " << _nX << " " << _nY << std::endl;
    for(int i=0; i<size(); ++i)
    {
        for(int j=0; j<_nX+_nY; ++j)
        {
            file << _data[i][j] << "\t";
        }
        file << std::endl;
    }
    file.close();
}

vec clustering::get(int i) const
{
    return _data[i];
}

vec clustering::operator[](int i) const
{
    return _data[i];
}

int clustering::size() const
{
    return _data.size();
}

vec clustering::min() const
{
    return _min;
}

vec clustering::max() const
{
    return _max;
}
