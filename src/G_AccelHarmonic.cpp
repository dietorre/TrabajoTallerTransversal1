//$Source$
//------------------------------------------------------------------------------
//      G_AccelHarmonic
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/11

#include "../include/G_AccelHarmonic.h"

/**
 * @brief Computes the gradient of the Earth's harmonic gravity field.
 *
 * @param r Satellite position vector in the true-of-date system.
 * @param U Transformation matrix to body-fixed system.
 * @param n Gravity model degree.
 * @param m Gravity model order.
 *
 * @return Gradient (G=da/dr) in the true-of-date system.
 */
Matrix G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max)
{
    double d = 1.0; // Position increment [m]

    Matrix G(3, 3);
    Matrix dr(3, 1);

    // Gradient
    for (int i = 1; i <= 3; i++)
    {
        // Set offset in i-th component of the position vector
        dr = Matrix(3, 1);
        dr(i) = d;
        // Acceleration difference
        Matrix da = AccelHarmonic(r + (1.0 / 2) * dr, U, n_max, m_max) - AccelHarmonic(r - (1.0 / 2) * dr, U, n_max, m_max);
        // Derivative with respect to i-th axis
        Matrix aux = (1 / d) * da;
        G(1, i) = aux(1);
        G(2, i) = aux(2);
        G(3, i) = aux(3);
    }

    return G;
}