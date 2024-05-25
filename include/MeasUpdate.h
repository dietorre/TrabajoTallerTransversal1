//$Header$
//------------------------------------------------------------------------------
//      MeasUpdate
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
 *   Updates the state estimate and covariance matrix using the measurement
 * **/

#ifndef MeasUpdate_h
#define MeasUpdate_h

#include "Matrix.h"

void MeasUpdate(Matrix &K, Matrix &x, double z, double g, double s, Matrix G, Matrix &P, int n);

#endif