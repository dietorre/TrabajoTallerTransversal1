//$Header$
//------------------------------------------------------------------------------
//      EccAnom
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/27


#include "../include/EccAnom.h"

/**
 * @brief Computes the eccentric anomaly for elliptic orbits.
 *
 * @param M Mean anomaly in [rad].
 * @param e Eccentricity of the orbit [0,1].
 *
 * @return Eccentric anomaly in [rad].
 *
 * @note Last modified:   2015/08/12   M. Mahooti
 */
double EccAnom(double M, double e)
{
    int maxit = 15;
    int i = 1;

    // Starting value
    M = fmod(M, Const::pi2);
    double E;

    if (e < 0.8)
    {
        E = M;
    }
    else
    {
        E = M_PI;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // double eps = std::numeric_limits<double>::epsilon();
    double eps = 2.22044604925031e-16;

    // Iteration
    while (fabs(f) > 1e2 * eps)
    {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit)
        {
            throw std::runtime_error("Problemas de convergencia en EccAnom");
        }
    }

    return E;
}