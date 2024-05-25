//$Header$
//------------------------------------------------------------------------------
//      NutAngles
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28
/**
 *
 * Purpose:
 *   Nutation in longitude and obliquity
 * **/

#ifndef NutAngles_h
#define NutAngles_h

#include "Matrix.h"
#include "SAT_Const.h"

void NutAngles(double &dpsi, double &deps, double Mjd_TT);

#endif