//$Source$
//------------------------------------------------------------------------------
//      NutMatrix
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/09

#include "../include/NutMatrix.h"

/**
 * @brief Computes the transformation from mean to true equator and equinox.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 *
 * @return NutMat Nutation matrix.
 */
Matrix NutMatrix(double Mjd_TT)
{
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles(dpsi, deps, Mjd_TT);

    // Transformation from mean to true equator and equinox
    return R_x(-eps - deps) * R_z(-dpsi) * R_x(+eps);
}
