//$Source$
//------------------------------------------------------------------------------
//      TimeUpdate
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/11

#include "../include/TimeUpdate.h"

/**
 * @brief Updates the covariance matrix using the time update step.
 *
 * @param P Covariance matrix before the update.
 * @param Phi State transition matrix.
 * @param Qdt Covariance matrix of the process noise.
 */
void TimeUpdate(Matrix &P, Matrix Phi, double Qdt)
{
    P = Phi*P*Phi.transponer() + Qdt;
}