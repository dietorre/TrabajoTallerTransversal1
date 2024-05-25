//$Source$
//------------------------------------------------------------------------------
//      AzElPa
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28

#include "../include/AzElPa.h"

/**
 * @brief Computes azimuth, elevation, and their partial derivatives from local tangent coordinates.
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame).
 *
 * @param[out] A Azimuth [rad].
 * @param[out] E Elevation [rad].
 * @param[out] dAds Partial derivatives of azimuth with respect to s.
 * @param[out] dEds Partial derivatives of elevation with respect to s.
 */
void AzElPa(double &Az, double &El, Matrix &dAds, Matrix &dEds, Matrix s)
{

    dAds = Matrix(1, 3);
    dEds = Matrix(1, 3);

    double pi2 = Const::pi2;
    double rho = sqrt(s(1) * s(1) + s(2) * s(2));

    // Angles
    Az = atan2(s(1), s(2));

    if (Az < 0.0)
    {
        Az = Az + pi2;
    }

    El = atan(s(3) / rho);

    // Partials
    dAds(1) = s(2) / (rho * rho);
    dAds(2) = -s(1) / (rho * rho);
    dAds(3) = 0.0;

    dEds(1) = -s(1) * s(3) / rho * (1 / dot(s.transponer(), s));
    dEds(2) = -s(2) * s(3) / rho * (1 / dot(s.transponer(), s));
    dEds(3) = rho * (1 / dot(s.transponer(), s));
}