//$Source$
//------------------------------------------------------------------------------
//      MeasUpdate
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/11

#include "../include/MeasUpdate.h"

/**
 * @brief Updates the state estimate and covariance matrix using the measurement.
 *
 * @param[out] K Kalman gain.
 * @param[in,out] x State vector.
 * @param[in] z Measurement.
 * @param[in] g Measurement model prediction.
 * @param[in] s Measurement noise standard deviation.
 * @param[in] G Measurement model Jacobian matrix.
 * @param[in,out] P Covariance matrix of the state estimate.
 * @param[in] n Number of state variables.
 *
 */
void MeasUpdate(Matrix & K, Matrix & x, double z, double g, double s, Matrix G, Matrix & P, int n)
{
    int m = 1;
    Matrix Inv_W(m,m);

    for(int i=1;i<=m;i++){
        Inv_W(i,i) = s*s;    // Inverse weight (measurement covariance)
    }

    // Kalman gain
    Matrix aux1 = P*G.transponer();
    Matrix aux2 = (Inv_W+G*P*G.transponer()).inverse();
    K = P*G.transponer()*(Inv_W+G*P*G.transponer()).inverse();

    // State update
    x = x + K*(z-g);

    // Covariance update
    Matrix identidad(n,n);
    for(int i=1;i<=n;i++){
        identidad(i,i) = 1;
    }
    P = (identidad-K*G)*P;
}