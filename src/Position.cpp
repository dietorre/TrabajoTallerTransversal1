//$Source$
//------------------------------------------------------------------------------
//      Position
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28

#include "../include/Position.h"

/**
 * @brief Computes the position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]).
 *
 * @param lon Geodetic East longitude [rad].
 * @param lat Geodetic latitude [rad].
 * @param h Altitude [m].
 *
 * @return r Position vector [m].
 *
 */
Matrix Position(double lon, double lat, double h)
{
    double R_equ = Const::R_Earth;
    double f = Const::f_Earth;

    double e2 = f * (2.0 - f); // Square of eccentricity
    double CosLat = cos(lat);  // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    Matrix r(3, 1);

    r(1) = (N + h) * CosLat * cos(lon);
    r(2) = (N + h) * CosLat * sin(lon);
    r(3) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}