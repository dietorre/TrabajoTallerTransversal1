//$Header$
//------------------------------------------------------------------------------
//      gast
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
 *   Greenwich Apparent Sidereal Time
 * **/

#ifndef gast_h
#define gast_h

#include "gmst.h"
#include "EqnEquinox.h"
#include <cmath>

double gast(double Mjd_UT1);

#endif