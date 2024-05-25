//$Header$
//------------------------------------------------------------------------------
//      GHAMatrix
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/10
/**
 *
 * Purpose:
 *   Transformation from true equator and equinox to Earth equator and 
 *   Greenwich meridian system 
 * **/

#ifndef GHA_Matrix
#define GHA_Matrix

#include "Matrix.h"
#include "R_z.h"
#include "gast.h"

Matrix GHAMatrix(double Mjd_UT1);


#endif