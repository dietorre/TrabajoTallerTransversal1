//$Header$
//------------------------------------------------------------------------------
//      VarEqn
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
 *   Computes the variational equations, i.e. the derivative of the state vector.
 *   and the state transition matrix
 * **/

#ifndef VarEqn_h
#define VarEqn_h

#include "Matrix.h"
#include "SAT_Const.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "parametros.h"

Matrix VarEqn(double x, Matrix yPhi);

#endif