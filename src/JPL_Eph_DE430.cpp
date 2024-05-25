//$Source$
//------------------------------------------------------------------------------
//      JPL_Eph_DE430
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/10

#include "../include/JPL_Eph_DE430.h"

/**
 * @brief Computes the sun, moon, and nine major planets' equatorial position using JPL Ephemerides (DE430).
 *
 * @param[out] r_Mercury Equatorial position of Mercury.
 * @param[out] r_Venus Equatorial position of Venus.
 * @param[out] r_Earth Equatorial position of the Earth (solar system barycenter (SSB)).
 * @param[out] r_Mars Equatorial position of Mars.
 * @param[out] r_Jupiter Equatorial position of Jupiter.
 * @param[out] r_Saturn Equatorial position of Saturn.
 * @param[out] r_Uranus Equatorial position of Uranus.
 * @param[out] r_Neptune Equatorial position of Neptune.
 * @param[out] r_Pluto Equatorial position of Pluto.
 * @param[out] r_Moon Equatorial position of the Moon.
 * @param[out] r_Sun Equatorial position of the Sun (geocentric equatorial position referred to the International Celestial Reference Frame (ICRF)).
 * @param[in] Mjd_TDB Modified Julian Date of TDB.
 *
 * @note Light-time is already taken into account.
 */
void JPL_Eph_DE430(Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun, double Mjd_TDB)
{
    extern Matrix PC;
    double JD = Mjd_TDB + 2400000.5;

    int index = -1;
    for (int i = 1; i <= PC.columnas(); ++i)
    {
        if (PC(i, 1) <= JD && JD <= PC(i, 2))
        {
            index = i;
            break;
        }
    }

    int i = index;

    Matrix PCtemp = PC.getFila(i);

    double t1 = PCtemp(1) - 2400000.5; // MJD at start of interval

    double dt = Mjd_TDB - t1;

    Matrix temp(231, 13, 270);

    Matrix Cx_Earth = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Earth = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Earth = PCtemp.slice(temp(3), temp(4) - 1);
    temp = temp + 39;

    Matrix Cx = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz = PCtemp.slice(temp(3), temp(4) - 1);

    Cx_Earth = Matrix(Cx_Earth, Cx);
    Cy_Earth = Matrix(Cy_Earth, Cy);
    Cz_Earth = Matrix(Cz_Earth, Cz);

    int j;
    double Mjd0;

    if (0 <= dt && dt <= 16)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (16 < dt && dt <= 32)
    {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    r_Earth = 1e3 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, Cx_Earth.slice(13 * j + 1, 13 * j + 13), Cy_Earth.slice(13 * j + 1, 13 * j + 13), Cz_Earth.slice(13 * j + 1, 13 * j + 13)).transponer();
    temp = Matrix(441, 13, 480);

    Matrix Cx_Moon = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Moon = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Moon = PCtemp.slice(temp(3), temp(4) - 1);

    for (int i = 1; i <= 7; i++)
    {
        temp = temp + 39;
        Cx = PCtemp.slice(temp(1), temp(2) - 1);
        Cy = PCtemp.slice(temp(2), temp(3) - 1);
        Cz = PCtemp.slice(temp(3), temp(4) - 1);
        Cx_Moon = Matrix(Cx_Moon, Cx);
        Cy_Moon = Matrix(Cy_Moon, Cy);
        Cz_Moon = Matrix(Cz_Moon, Cz);
    }

    if (0 <= dt && dt <= 4)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (4 < dt && dt <= 8)
    {
        j = 1;
        Mjd0 = t1 + 4 * j;
    }
    else if (8 < dt && dt <= 12)
    {
        j = 2;
        Mjd0 = t1 + 4 * j;
    }
    else if (12 < dt && dt <= 16)
    {
        j = 3;
        Mjd0 = t1 + 4 * j;
    }
    else if (16 < dt && dt <= 20)
    {
        j = 4;
        Mjd0 = t1 + 4 * j;
    }

    else if (20 < dt && dt <= 24)
    {
        j = 5;
        Mjd0 = t1 + 4 * j;
    }

    else if (24 < dt && dt <= 28)
    {
        j = 6;
        Mjd0 = t1 + 4 * j;
    }

    else if (28 < dt && dt <= 32)
    {
        j = 7;
        Mjd0 = t1 + 4 * j;
    }

    r_Moon = 1e3 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, Cx_Moon.slice(13 * j + 1, 13 * j + 13), Cy_Moon.slice(13 * j + 1, 13 * j + 13), Cz_Moon.slice(13 * j + 1, 13 * j + 13)).transponer();

    temp = Matrix(753, 11, 786);

    Matrix Cx_Sun = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Sun = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Sun = PCtemp.slice(temp(3), temp(4) - 1);
    temp = temp + 33;
    Cx = PCtemp.slice(temp(1), temp(2) - 1);
    Cy = PCtemp.slice(temp(2), temp(3) - 1);
    Cz = PCtemp.slice(temp(3), temp(4) - 1);
    Cx_Sun = Matrix(Cx_Sun, Cx);
    Cy_Sun = Matrix(Cy_Sun, Cy);
    Cz_Sun = Matrix(Cz_Sun, Cz);
    if (0 <= dt && dt <= 16)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (16 < dt && dt <= 32)
    {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    r_Sun = 1e3 * Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, Cx_Sun.slice(11 * j + 1, 11 * j + 11), Cy_Sun.slice(11 * j + 1, 11 * j + 11), Cz_Sun.slice(11 * j + 1, 11 * j + 11)).transponer();

    temp = Matrix(3, 14, 45);
    Matrix Cx_Mercury = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Mercury = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Mercury = PCtemp.slice(temp(3), temp(4) - 1);
    for (int i = 1; i <= 3; i++)
    {
        temp = temp + 42;
        Cx = PCtemp.slice(temp(1), temp(2) - 1);
        Cy = PCtemp.slice(temp(2), temp(3) - 1);
        Cz = PCtemp.slice(temp(3), temp(4) - 1);
        Cx_Mercury = Matrix(Cx_Mercury, Cx);
        Cy_Mercury = Matrix(Cy_Mercury, Cy);
        Cz_Mercury = Matrix(Cz_Mercury, Cz);
    }
    if (0 <= dt && dt <= 8)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (8 < dt && dt <= 16)
    {
        j = 1;
        Mjd0 = t1 + 8 * j;
    }
    else if (16 < dt && dt <= 24)
    {
        j = 2;
        Mjd0 = t1 + 8 * j;
    }
    else if (24 < dt && dt <= 32)
    {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    r_Mercury = 1e3 * Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, Cx_Mercury.slice(14 * j + 1, 14 * j + 14), Cy_Mercury.slice(14 * j + 1, 14 * j + 14), Cz_Mercury.slice(14 * j + 1, 14 * j + 14)).transponer();

    temp = Matrix(171, 10, 201);
    Matrix Cx_Venus = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Venus = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Venus = PCtemp.slice(temp(3), temp(4) - 1);
    temp = temp + 30;
    Cx = PCtemp.slice(temp(1), temp(2) - 1);
    Cy = PCtemp.slice(temp(2), temp(3) - 1);
    Cz = PCtemp.slice(temp(3), temp(4) - 1);
    Cx_Venus = Matrix(Cx_Venus, Cx);
    Cy_Venus = Matrix(Cy_Venus, Cy);
    Cz_Venus = Matrix(Cz_Venus, Cz);
    if (0 <= dt && dt <= 16)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (16 < dt && dt <= 32)
    {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    r_Venus = 1e3 * Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16, Cx_Venus.slice(10 * j + 1, 10 * j + 10), Cy_Venus.slice(10 * j + 1, 10 * j + 10), Cz_Venus.slice(10 * j + 1, 10 * j + 10)).transponer();

    temp = Matrix(309, 11, 342);
    Matrix Cx_Mars = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Mars = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Mars = PCtemp.slice(temp(3), temp(4) - 1);
    j = 0;
    Mjd0 = t1;
    r_Mars = 1e3 * Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, Cx_Mars.slice(11 * j + 1, 11 * j + 11), Cy_Mars.slice(11 * j + 1, 11 * j + 11), Cz_Mars.slice(11 * j + 1, 11 * j + 11)).transponer();

    temp = Matrix(342, 8, 366);
    Matrix Cx_Jupiter = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Jupiter = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Jupiter = PCtemp.slice(temp(3), temp(4) - 1);
    j = 0;
    Mjd0 = t1;
    r_Jupiter = 1e3 * Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, Cx_Jupiter.slice(8 * j + 1, 8 * j + 8), Cy_Jupiter.slice(8 * j + 1, 8 * j + 8), Cz_Jupiter.slice(8 * j + 1, 8 * j + 8)).transponer();

    temp = Matrix(366, 7, 387);
    Matrix Cx_Saturn = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Saturn = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Saturn = PCtemp.slice(temp(3), temp(4) - 1);
    j = 0;
    Mjd0 = t1;
    r_Saturn = 1e3 * Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, Cx_Saturn.slice(7 * j + 1, 7 * j + 7), Cy_Saturn.slice(7 * j + 1, 7 * j + 7), Cz_Saturn.slice(7 * j + 1, 7 * j + 7)).transponer();

    temp = Matrix(387, 6, 405);
    Matrix Cx_Uranus = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Uranus = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Uranus = PCtemp.slice(temp(3), temp(4) - 1);
    j = 0;
    Mjd0 = t1;
    r_Uranus = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Uranus.slice(6 * j + 1, 6 * j + 6), Cy_Uranus.slice(6 * j + 1, 6 * j + 6), Cz_Uranus.slice(6 * j + 1, 6 * j + 6)).transponer();

    temp = Matrix(405, 6, 423);
    Matrix Cx_Neptune = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Neptune = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Neptune = PCtemp.slice(temp(3), temp(4) - 1);
    j = 0;
    Mjd0 = t1;
    r_Neptune = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Neptune.slice(6 * j + 1, 6 * j + 6), Cy_Neptune.slice(6 * j + 1, 6 * j + 6), Cz_Neptune.slice(6 * j + 1, 6 * j + 6)).transponer();

    temp = Matrix(423, 6, 441);
    Matrix Cx_Pluto = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Pluto = PCtemp.slice(temp(2), temp(3) - 1);
    Matrix Cz_Pluto = PCtemp.slice(temp(3), temp(4) - 1);
    j = 0;
    Mjd0 = t1;
    r_Pluto = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Pluto.slice(6 * j + 1, 6 * j + 6), Cy_Pluto.slice(6 * j + 1, 6 * j + 6), Cz_Pluto.slice(6 * j + 1, 6 * j + 6)).transponer();

    temp = Matrix(819, 10, 839);
    Matrix Cx_Nutations = PCtemp.slice(temp(1), temp(2) - 1);
    Matrix Cy_Nutations = PCtemp.slice(temp(2), temp(3) - 1);
    for (int i = 1; i <= 3; i++)
    {
        temp = temp + 20;
        Cx = PCtemp.slice(temp(1), temp(2) - 1);
        Cy = PCtemp.slice(temp(2), temp(3) - 1);
        Cx_Nutations = Matrix(Cx_Nutations, Cx);
        Cy_Nutations = Matrix(Cy_Nutations, Cy);
    }
    if (0 <= dt && dt <= 8)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (8 < dt && dt <= 16)
    {
        j = 1;
        Mjd0 = t1 + 8 * j;
    }
    else if (16 < dt && dt <= 24)
    {
        j = 2;
        Mjd0 = t1 + 8 * j;
    }
    else if (24 < dt && dt <= 32)
    {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    // Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations(10*j+1,10*j+10),Cy_Nutations(10*j+1,10*j+10),zeros(10,1)).transponer();

    // temp = (899,10,929);
    // Cx_Librations = PCtemp(temp(1),temp(2)-1);
    // Cy_Librations = PCtemp(temp(2),temp(3)-1);
    // Cz_Librations = PCtemp(temp(3),temp(4)-1);
    // for i=1,3
    //     temp = temp+30;
    //     Cx = PCtemp(temp(1),temp(2)-1);
    //     Cy = PCtemp(temp(2),temp(3)-1);
    //     Cz = PCtemp(temp(3),temp(4)-1);
    //     Cx_Librations = [Cx_Librations,Cx];
    //     Cy_Librations = [Cy_Librations,Cy];
    //     Cz_Librations = [Cz_Librations,Cz];
    // end
    // if (0<=dt && dt<=8)
    //     j=0;
    //     Mjd0 = t1;
    // else if(8<dt && dt<=16)
    //     j=1;
    //     Mjd0 = t1+8*j;
    // else if (16<dt && dt<=24)
    //     j=2;
    //     Mjd0 = t1+8*j;
    // else if(24<dt && dt<=32)
    //     j=3;
    //     Mjd0 = t1+8*j;
    // end
    // Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations(10*j+1,10*j+10),...
    //                     Cy_Librations(10*j+1,10*j+10), Cz_Librations(10*j+1,10*j+10)).transponer();
    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1 / (1 + EMRAT);
    r_Earth = r_Earth - EMRAT1 * r_Moon;
    r_Mercury = -1 * r_Earth + r_Mercury;
    r_Venus = -1 * r_Earth + r_Venus;
    r_Mars = -1 * r_Earth + r_Mars;
    r_Jupiter = -1 * r_Earth + r_Jupiter;
    r_Saturn = -1 * r_Earth + r_Saturn;
    r_Uranus = -1 * r_Earth + r_Uranus;
    r_Neptune = -1 * r_Earth + r_Neptune;
    r_Pluto = -1 * r_Earth + r_Pluto;
    r_Sun = -1 * r_Earth + r_Sun;
}