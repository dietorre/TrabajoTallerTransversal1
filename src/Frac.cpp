//$Source$
//------------------------------------------------------------------------------
//      Frac
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/27

#include "../include/Frac.h"
#include <cmath>

/**
 * @brief Fractional part of a number (y=x-[x]).
 *
 * @param x Input number.
 *
 * @return Fractional part of the input number.
 *
 */
double Frac(double x)
{
    return x - floor(x);
}