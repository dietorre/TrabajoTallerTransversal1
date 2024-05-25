//$Header$
//------------------------------------------------------------------------------
//      Position
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
 *   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
 *   latitude [rad], altitude [m])
 * **/


#ifndef Position_h
#define Position_h

#include "Matrix.h"
#include "SAT_Const.h"

Matrix Position(double lon, double lat, double h);

#endif