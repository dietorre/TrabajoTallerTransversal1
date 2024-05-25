//$Header$
//------------------------------------------------------------------------------
//      PoleMatrix
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/09
/**
 *
 * Purpose:
 *   Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 * **/

#ifndef Pole_Matrix
#define Pole_Matrix

#include "Matrix.h"
#include "R_x.h"
#include "R_y.h"

Matrix PoleMatrix(double xp, double yp);

#endif