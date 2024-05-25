//$Header$
//------------------------------------------------------------------------------
//      PrecMatrix
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
 *   Precession transformation of equatorial coordinates
 * **/

#ifndef Prec_Matrix
#define Prec_Matrix

#include "Matrix.h"
#include "SAT_Const.h"
#include "R_z.h"
#include "R_x.h"
#include "R_y.h"

Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif