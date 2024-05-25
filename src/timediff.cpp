//$Header$
//------------------------------------------------------------------------------
//      timediff
//------------------------------------------------------------------------------
// Proyecto TT1
//
// **Legal**
//
// Author: Diego de la Torre
// Created: 2024/04/24

#include "../include/timediff.h"

/**
 * @brief Computes various time differences.
 *
 * @param UT1_UTC   UT1 - UTC [s].
 * @param TAI_UTC   TAI - UTC [s].
 * @param[out] UT1_TAI UT1 - TAI [s].
 * @param[out] UTC_GPS UTC - GPS [s].
 * @param[out] UT1_GPS UT1 - GPS [s].
 * @param[out] TT_UTC  TT - UTC [s].
 * @param[out] GPS_UTC GPS - UTC [s].
 */
void timediff(double UT1_UTC, double TAI_UTC, double &UT1_TAI, double &UTC_GPS, double &UT1_GPS, double &TT_UTC, double &GPS_UTC)
{
    double TT_TAI = +32.184; // TT-TAI time difference [s]

    double GPS_TAI = -19.0; // GPS-TAI time difference [s]

    double TT_GPS = TT_TAI - GPS_TAI; // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI; // TAI-GPS time difference [s]

    UT1_TAI = UT1_UTC - TAI_UTC; // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC; // UTC-TAI time difference [s]

    UTC_GPS = UTC_TAI - GPS_TAI; // UTC_GPS time difference [s]

    UT1_GPS = UT1_TAI - GPS_TAI; // UT1-GPS time difference [s]

    TT_UTC = TT_TAI - UTC_TAI; //  TT-UTC time difference [s]

    GPS_UTC = GPS_TAI - UTC_TAI; // GPS-UTC time difference [s]
}