//$Header$
//------------------------------------------------------------------------------
//      AzElPa
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
 *   Computes azimuth, elevation and partials from local tangent coordinates
 * **/

#ifndef AzElPa_h
#define AzElPa_h

#include "Matrix.h"
#include "SAT_Const.h"

void AzElPa(double &Az, double &El, Matrix &dAds, Matrix &dEds, Matrix s);

#endif