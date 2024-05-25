//$Source$
//------------------------------------------------------------------------------
//      Geodetic
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28

#include "../include/Geodetic.h"

/**
 * @brief Computes geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 * from a given position vector (r [m]).
 *
 * @param r Position vector [m].
 *
 * @return lon Longitude [rad].
 * @return lat Latitude [rad].
 * @return h Altitude [m].
 *
 */

void Geodetic(double &lon, double &lat, double &h, Matrix r)
{
    double R_equ = Const::R_Earth;
    double f = Const::f_Earth;

    // double eps = std::numeric_limits<double>::epsilon();
    double eps = 2.22044604925031e-16;

    double epsRequ = eps * R_equ; // Convergence criterion
    double e2 = f * (2.0 - f);    // Square of eccentricity

    double X = r(1); // Cartesian coordinates
    double Y = r(2);
    double Z = r(3);
    double rho2 = X * X + Y * Y; // Square of distance from z-axis

    // Check validity of input data
    if (norm(r) == 0.0)
    {
        throw std::invalid_argument("invalid input in Geodetic constructor");
        lon = 0.0;
        lat = 0.0;
        h = -R_equ;
    }

    // Iteration
    double dZ = e2 * Z;
    double ZdZ = 0;
    double Nh = 0;
    double N = 0;

    while (1)
    {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        double SinPhi = ZdZ / Nh; // Sine of geodetic latitude
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        double dZ_new = N * e2 * SinPhi;
        if (fabs(dZ - dZ_new) < epsRequ)
        {
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2(Y, X);
    lat = atan2(ZdZ, sqrt(rho2));
    h = Nh - N;
}