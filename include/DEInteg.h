//$Header$
//------------------------------------------------------------------------------
//      DEInteg
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
 *   Numerical integration methods for ordinaray differential equations
 *
 *   This module provides implemenation of the variable order variable
 *   stepsize multistep method of Shampine & Gordon.
 * **/

#ifndef DEInteg_h
#define DEInteg_h

#include "Matrix.h"
#include "sign_.h"
#include <cmath>

void DEInteg(Matrix (*func)(double, Matrix), double t, double tout, double relerr, double abserr, double n_eqn, Matrix &y);

#endif