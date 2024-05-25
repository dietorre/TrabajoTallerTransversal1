//$Source$
//------------------------------------------------------------------------------
//      GHAMatrix
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/10

#include "../include/GHAMatrix.h"

/**
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system.
 *
 * @param Mjd_UT1 Modified Julian Date UT1.
 *
 * @return GHAmat Greenwich Hour Angle matrix.
 */

Matrix GHAMatrix(double Mjd_UT1)
{
    return R_z(gast(Mjd_UT1));
}