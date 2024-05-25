//$Header$
//------------------------------------------------------------------------------
//      TimeUpdate
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
 *   Updates the covariance matrix using the time update step
 * **/

#ifndef TimeUpdate_h
#define TimeUpdate_h

#include "Matrix.h"

void TimeUpdate(Matrix &P, Matrix Phi, double Qdt = 0);

#endif