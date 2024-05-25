//$Source$
//------------------------------------------------------------------------------
//      mjday
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/21

#include "../include/Mjday.h"
#include <cmath>

/**
 * @brief Computes the Modified Julian Date.
 *
 * @param yr Year.
 * @param mon Month.
 * @param day Day.
 * @param hr Universal time hour.
 * @param min Universal time minute.
 * @param sec Universal time second.
 *
 * @return Mjd Modified Julian Date.
 */
double mjday(int yr, int mon, int day, int hr, int min, double sec)
{
	double jd = 367.0 * yr - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25) + floor(275 * mon / 9.0) + day + 1721013.5 + ((sec / 60.0 + min) / 60.0 + hr) / 24.0;

	return jd - 2400000.5;
}
