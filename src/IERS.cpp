//$ource$
//------------------------------------------------------------------------------
//      IERS
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/27

#include "../include/IERS.h"
#include <iostream>

using namespace std;

/**
 * @brief Management of IERS time and polar motion data
 *
 * @param x_pole CIP x-coordinate [arcsec].
 * @param y_pole CIP y-coordinate [arcsec].
 * @param UT1_UTC UT1-UTC [s].
 * @param LOD Length of day [s].
 * @param dpsi Celestial pole offset (Δψ) [arcsec].
 * @param deps Celestial pole offset (Δε) [arcsec].
 * @param dx_pole CIP x-coordinate rate [arcsec/s].
 * @param dy_pole CIP y-coordinate rate [arcsec/s].
 * @param TAI_UTC TAI-UTC [s].
 * @param eop Earth orientation parameters.
 * @param Mjd_UTC Modified Julian Date UTC.
 * @param interp Interpolation flag (default 'n').
 */
void IERS(double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC, Matrix eop, double Mjd_UTC, char interp)
{

    if (interp == 'l')
    {
        // linear interpolation
        float mjd = (floor(Mjd_UTC));
        int i = 0;

        for (i = 1; i <= eop.columnas(); i++)
        {
            if (mjd == eop(4, i))
            {
                break;
            }
        }

        Matrix preeop = eop.getColumna(i);
        Matrix nexteop = eop.getColumna(i + 1);
        double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = preeop(5) + (nexteop(5) - preeop(5)) * fixf;
        y_pole = preeop(6) + (nexteop(6) - preeop(6)) * fixf;
        UT1_UTC = preeop(7) + (nexteop(7) - preeop(7)) * fixf;
        LOD = preeop(8) + (nexteop(8) - preeop(8)) * fixf;
        dpsi = preeop(9) + (nexteop(9) - preeop(9)) * fixf;
        deps = preeop(10) + (nexteop(10) - preeop(10)) * fixf;
        dx_pole = preeop(11) + (nexteop(11) - preeop(11)) * fixf;
        dy_pole = preeop(12) + (nexteop(12) - preeop(12)) * fixf;
        TAI_UTC = preeop(13);

        x_pole = x_pole / Const::Arcs; // Pole coordinate [rad]
        y_pole = y_pole / Const::Arcs; // Pole coordinate [rad]
        dpsi = dpsi / Const::Arcs;
        deps = deps / Const::Arcs;
        dx_pole = dx_pole / Const::Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole / Const::Arcs; // Pole coordinate [rad]
    }

    else if (interp == 'n')
    {

        float mjd = (floor(Mjd_UTC));

        int i = 1;

        for (i = 1; i <= eop.columnas(); i++)
        {
            if (mjd == eop(4, i))
            {
                break;
            }
        }

        Matrix preeop = eop.getColumna(i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = preeop(5) / Const::Arcs; // Pole coordinate [rad]
        y_pole = preeop(6) / Const::Arcs; // Pole coordinate [rad]
        UT1_UTC = preeop(7);              // UT1-UTC time difference [s]
        LOD = preeop(8);                  // Length of day [s]
        dpsi = preeop(9) / Const::Arcs;
        deps = preeop(10) / Const::Arcs;
        dx_pole = preeop(11) / Const::Arcs; // Pole coordinate [rad]
        dy_pole = preeop(12) / Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = preeop(13);               // TAI-UTC time difference [s]
    }
}