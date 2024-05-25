//$Header$
//------------------------------------------------------------------------------
//      Geodetic
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
 *   Greenwich Apparent Sidereal Time
 * **/

#ifndef Geodetic_H
#define Geodetic_H

#include "Matrix.h"
#include "SAT_Const.h"

void Geodetic(double &lon, double &lat, double &h, Matrix r);

#endif 