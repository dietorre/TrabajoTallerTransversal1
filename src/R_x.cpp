//$Header$
//------------------------------------------------------------------------------
//      R_x
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/11

#include "../include/R_x.h"
#include "../include/Matrix.h"
#include <cmath>

/**
 * @brief Computes the rotation matrix around the x-axis.
 *
 * @param angle Angle of rotation [rad].
 *
 * @return rotmat Rotation matrix.
 */
Matrix R_x(double angle)
{
    double C, S;
    Matrix rotmat(3, 3);

    C = cos(angle);
    S = sin(angle);

    rotmat(1, 1) = 1.0;
    rotmat(1, 2) = 0.0;
    rotmat(1, 3) = 0.0;
    rotmat(2, 1) = 0.0;
    rotmat(2, 2) = C;
    rotmat(2, 3) = S;
    rotmat(3, 1) = 0.0;
    rotmat(3, 2) = -1.0 * S;
    rotmat(3, 3) = C;

    return rotmat;
}