//$Source$
//------------------------------------------------------------------------------
//      R_z
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/11


#include "../include/R_z.h"
#include "../include/Matrix.h"
#include <cmath>

/**
 * @brief Computes the rotation matrix around the z-axis
 *
 * @param angle Angle of rotation [rad].
 *
 * @return rotmat Resulting rotation vector.
 */
Matrix R_z(double angle){
    double C, S;
    Matrix rotmat(3,3);

    C = cos(angle);
    S = sin(angle);

    rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
    rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
    rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;

    return rotmat;
}