//$Source$
//------------------------------------------------------------------------------
//      MeanObliquity
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28

#include "../include/MeanObliquity.h"

/**
 * @brief Computes the mean obliquity of the ecliptic.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 *
 * @return MOblq Mean obliquity of the ecliptic [rad].
 */
double MeanObliquity(double Mjd_TT)
{
    double T = (Mjd_TT - Const::MJD_J2000) / 36525;

    return Const::Rad * (84381.448 / 3600 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600);
}