//$Source$
//------------------------------------------------------------------------------
//      PoleMatrix
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/09

#include "../include/PoleMatrix.h"

/**
 * @brief Computes the transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 * @param xp Pole coordinate along the x-axis.
 * @param yp Pole coordinate along the y-axis.
 *
 * @return PoleMat Pole matrix.
 */
Matrix PoleMatrix(double xp, double yp)
{
    return R_y(-xp) * R_x(-yp);
}