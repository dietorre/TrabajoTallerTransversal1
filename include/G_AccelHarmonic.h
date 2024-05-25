//$Header$
//------------------------------------------------------------------------------
//      G_AccelHarmonic
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
 *   Computes the gradient of the Earth's harmonic gravity field
 * **/

#ifndef G_AccelHarmonic_h
#define G_AccelHarmonic_h

#include "Matrix.h"
#include "AccelHarmonic.h"

Matrix G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max);

#endif