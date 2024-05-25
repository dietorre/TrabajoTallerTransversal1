//$Source$
//------------------------------------------------------------------------------
//      gast
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/10

#include "../include/gast.h"

/**
 * @brief Computes Greenwich Apparent Sidereal Time.
 *
 * @param Mjd_UT1 Modified Julian Date UT1.
 *
 * @return GAST in [rad].
 */
double gast(double Mjd_UT1)
{
    return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * M_PI);
}