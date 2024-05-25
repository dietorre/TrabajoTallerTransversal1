//$Source$
//------------------------------------------------------------------------------
//      sign_
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/24

#include "../include/sign_.h"
#include "cmath"

/**
 * @brief Returns the absolute value of a with the sign of b.
 *
 * @param a Input value.
 * @param b Reference value.
 *
 * @return result Absolute value of a with the sign of b.
 */
double sign_(double a, double b)
{
    if (b >= 0.0)
        return fabs(a);
    else
        return -fabs(a);
}
