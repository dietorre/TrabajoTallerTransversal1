//$Header$
//------------------------------------------------------------------------------
//      JPL_Eph_DE430
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/05/10
/**
 *
 * Purpose:
 *   Computes the sun, moon, and nine major planets' equatorial
 *   position using JPL Ephemerides
 * **/

#ifndef JPL_Eph_DE430_h
#define JPL_Eph_DE430_h

#include "Matrix.h"
#include "Cheb3D.h"

void JPL_Eph_DE430(Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun, double Mjd_TDB);

#endif