//$Header$
//------------------------------------------------------------------------------
//      LTC
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/11
/**
 *
 * Purpose:
 *   Transformation from Greenwich meridian system to local tangent coordinates
 * **/

#ifndef LTC_h
#define LTC_h

#include "Matrix.h"
#include "R_y.h"
#include "R_z.h"

Matrix LTC(double lon, double lat);

#endif