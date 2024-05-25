//$Header$
//------------------------------------------------------------------------------
//      Cheb3D
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
 *   Chebyshev approximation of 3-dimensional vectors
 * **/

#ifndef Chab3D_h
#define Chab3D_h

#include "Matrix.h"
#include <iomanip>

Matrix Cheb3D(double t, double N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz);

#endif