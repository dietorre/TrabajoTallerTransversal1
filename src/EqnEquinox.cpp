//$Source$
//------------------------------------------------------------------------------
//      EqnEquinox
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/10

#include "../include/EqnEquinox.h"

/**
 * @brief Computes the equation of the equinoxes.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 *
 * @return Equation of the equinoxes.
 *
 * @note
 * The equation of the equinoxes dpsi*cos(eps) is the right ascension of
 * the mean equinox referred to the true equator and equinox and is equal
 * to the difference between apparent and mean sidereal time.
 */
double EqnEquinox(double Mjd_TT)
{
    double dpsi, deps;
    // Nutation in longitude and obliquity
    NutAngles(dpsi, deps, Mjd_TT);

    // Equation of the equinoxes
    return dpsi * cos(MeanObliquity(Mjd_TT));
}