#pragma once

#include "data.h"

/*! \brief a static class allowing to change from one parametrization
 *  to another.
 *  \ingroup core
 *
 *  \details
 *  Any \typedef function object or \typedef data object should have
 *  an associated parametrization.
 */
class params
{
    //! \brief list of all supported parametrization. An unsupported
    //! parametrization will go under the name <em>unknown</em>.
    enum type
    {
        unknown
    };

    //! \brief static function for input type convertion.
    static vec convert(const vec& invec, params::type intype, params::type outtype)
    {
        return invec;
    }
};
