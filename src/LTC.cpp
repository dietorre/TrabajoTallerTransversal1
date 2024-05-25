//$Source$
//------------------------------------------------------------------------------
//      LTC
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/11

#include "../include/LTC.h"

/**
 * @brief Transformation from Greenwich meridian system to local tangent coordinates.
 *
 * @param[in] lon Geodetic East longitude [rad].
 * @param[in] lat Geodetic latitude [rad].
 *
 * @return M Rotation matrix from the Earth equator and Greenwich meridian
 *           to the local tangent (East-North-Zenith) coordinate system.
 */
Matrix LTC(double lon, double lat)
{
    Matrix M = R_y(-1.0 * lat) * R_z(lon);

    for (int j = 1; j <= 3; j++)
    {
        double Aux = M(1, j);
        M(1, j) = M(2, j);
        M(2, j) = M(3, j);
        M(3, j) = Aux;
    }

    return M;
}