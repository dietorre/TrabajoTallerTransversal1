//$Header$
//------------------------------------------------------------------------------
//      unit
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/28

#include "../include/unit.h"

/**
 * @brief Calculates a unit vector given the original vector. If a zero vector is input, the output vector is also set to zero.
 *
 * @param vec Input vector.
 *
 * @return outvec Unit vector.
 */
Matrix unit(Matrix vec)
{
    double small = 0.000001;
    double magv = norm(vec);
    Matrix outvec(vec.filas(), vec.columnas());

    if (magv > small)
        for (int i = 1; i <= 3; i++)
        {
            outvec(i) = vec(i) / magv;
        }
    else
    {
        for (int i = 1; i <= 3; i++)
        {
            outvec(i) = 0.0;
        }
    }

    return outvec;
}