//$Header$
//------------------------------------------------------------------------------
//      NutMatrix
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
 *   Transformation from mean to true equator and equinox
 * **/

#ifndef Nut_Matrix
#define Nut_Matrix

#include "Matrix.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include "R_z.h"
#include "R_x.h"
#include "R_y.h"

Matrix NutMatrix(double Mjd_TT);

#endif
