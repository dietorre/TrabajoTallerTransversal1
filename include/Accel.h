//$Header$
//------------------------------------------------------------------------------
//      Accel
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
 *   Computes the acceleration of an Earth orbiting satellite due to 
 *    - the Earth's harmonic gravity field, 
 *    - the gravitational perturbations of the Sun and Moon
 *    - the solar radiation pressure and
 *    - the atmospheric drag
 * 
 * **/

#ifndef Accel_h
#define Accel_h

#include "Matrix.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "SAT_Const.h"
#include "parametros.h"


Matrix Accel(double x, Matrix Y);

#endif