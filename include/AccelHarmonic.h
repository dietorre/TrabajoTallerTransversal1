//$Header$
//------------------------------------------------------------------------------
//      AccelHarmonic
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
 *   Computes the acceleration due to the harmonic gravity field of the 
 *   central body 
 * 
 * **/
#ifndef Accel_Harmonic
#define Accel_Harmonic

#include "Matrix.h"
#include "SAT_Const.h"
#include "Legendre.h"

Matrix AccelHarmonic(Matrix r, Matrix E, int n_max, int m_max);

#endif