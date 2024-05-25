//$Header$
//------------------------------------------------------------------------------
//      Cheb3D
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28

#include "../include/Cheb3D.h"

/**
 * @brief Chebyshev approximation of 3-dimensional vectors.
 *
 * @param N Number of coefficients.
 * @param Ta Begin interval.
 * @param Tb End interval.
 * @param Cx Coefficients of Chebyshev polynomial (x-coordinate).
 * @param Cy Coefficients of Chebyshev polynomial (y-coordinate).
 * @param Cz Coefficients of Chebyshev polynomial (z-coordinate).
 */
Matrix Cheb3D(double t, double N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz)
{
    // Check validity
    if ((t < Ta) || (Tb < t))
    {
        throw std::invalid_argument("ERROR: Time out of range in Cheb3D::Value");
    }
    // Clenshaw algorithm
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1(1, 3);
    Matrix f2(1, 3);

    for (int i = N; i >= 2; i--)
    {
        Matrix old_f1 = f1;
        f1(1) = Cx(i);
        f1(2) = Cy(i);
        f1(3) = Cz(i);
        f1 = 2 * tau * old_f1 - f2 + f1;
        f2 = old_f1;
    }

    Matrix c1(1, 3);
    c1(1) = Cx(1);
    c1(2) = Cy(1);
    c1(3) = Cz(1);

    Matrix aux = f1 * tau - f2 + c1;
    return aux;
}