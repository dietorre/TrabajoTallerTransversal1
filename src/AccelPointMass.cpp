//$Source$
//------------------------------------------------------------------------------
//      AccelPointMass
//------------------------------------------------------------------------------
// Proyecto TT1
// 
// **Legal** 
// 
// Author: Diego de la Torre
// Created: 2024/04/27

#include "../include/AccelPointMass.h"

/**
 * @brief Computes the perturbational acceleration due to a point mass.
 *
 * @param r Satellite position vector.
 * @param s Point mass position vector.
 * @param GM Gravitational coefficient of the point mass.
 * 
 * @return Acceleration (a=d^2r/dt^2).
 */
Matrix AccelPointMass(Matrix r, Matrix s, double GM)
{
    
    // Relative position vector of satellite w.r.t. point mass 
    Matrix d = r - s;

    // Acceleration 

    return ( d*(1/pow(d.norm(),3)) +  s *(1/pow(s.norm(),3)) ) * (-GM);

}