//$Header$
//------------------------------------------------------------------------------
//      IERS
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
 *   Management of IERS time and polar motion data
 * **/

#ifndef IERS_H
#define IERS_H

#include "Matrix.h"
#include "SAT_Const.h"
#include <cmath>

void IERS(double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC, Matrix eop, double Mjd_UTC,char interp = 'n');

#endif 