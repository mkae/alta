//#include "rational_function.h"

#include <core/common.h>
#include <core/rational_function.h>

ALTA_DLL_EXPORT function* provide_function()
{
    return new rational_function();
}

