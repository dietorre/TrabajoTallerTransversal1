//$Header$
//------------------------------------------------------------------------------
//      AccelPointMass
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/27
/**
 *
 * Purpose:
 *   Computes the perturbational acceleration due to a point
 * **/

#ifndef AccelPointMass_H
#define AccelPointMass_H

#include "Matrix.h"
#include <cmath>

Matrix AccelPointMass(Matrix r, Matrix s, double GM);

#endif