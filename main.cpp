#include <iostream>
#include "include/Matrix.h"
#include "include/Mjday.h"
#include "include/globales.h"
#include "include/initGlobales.h"
#include "include/DEInteg.h"
#include "include/Accel.h"
#include "include/LTC.h"
#include "include/SAT_Const.h"
#include "include/VarEqn.h"
#include "include/Position.h"
#include "include/TimeUpdate.h"
#include "include/AzElPa.h"
#include "include/MeasUpdate.h"

using namespace std;

int main(int, char**){

    // declaro las variables externas
    extern Matrix eopdata;
    extern Matrix Cnm;
    extern Matrix Snm;
    extern Matrix PC;
    extern Matrix obs;

    extern Parametros AuxParam;

    //las inicializo
    eop19620101();
    GGM03S();
    DE430Coeff();
    GEOS3();

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Const::Rad; // [rad]
    double sigma_el = 0.0139*Const::Rad; // [rad]

    // Kaena Point station
    double lat = Const::Rad*21.5748;     // [rad]
    double lon = Const::Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt);

    Matrix r2(3,1);
    r2(1) = 6221397.628578685224056243896484375;
    r2(2) = 2867713.77965737879276275634765625;
    r2(3) = 3006155.98509948886930942535400390625;

    Matrix v2(3,1);
    v2(1) = 4645.0472516180598177015781402587890625;
    v2(2) = -2752.2159158820422817370854318141937255859375;
    v2(3) = -7507.99940987030640826560556888580322265625;

    Matrix Y0_apr(6,1);
    Y0_apr(1) = r2(1);
    Y0_apr(2) = r2(2);
    Y0_apr(3) = r2(3);
    Y0_apr(4) = v2(1);
    Y0_apr(5) = v2(2);
    Y0_apr(6) = v2(3);

 

    double Mjd0 = mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    int n_eqn  = 6;

    DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);
    Matrix Y = Y0_apr;

    Matrix P(6,6);

    for(int i=1;i<=3;i++){
        P(i,i)=1e8;
    }

    for(int i=4;i<=6;i++){
        P(i,i)=1e3;
    }

    Matrix LT = LTC(lon,lat);

    Matrix yPhi(42,1);
    Matrix Phi(6,6);

    // Measurement loop
    double t = 0;

    for (int i = 1; i <= nobs; i++)
    {
        // cout << i << endl;
        // Previous step
        double t_old = t;
        Matrix Y_old = Y;
        
        // Time increment and propagation
        Mjd_UTC = obs(i,1);                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]

        double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;     
        IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');
        
        double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
        timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
        double Mjd_TT = Mjd_UTC + TT_UTC/86400;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        for(int ii=1;ii<=6;ii++){
            yPhi(ii) = Y_old(ii);
            for(int j=1;j<=6;j++){
                if (ii==j){ 
                    yPhi(6*j+ii) = 1; 
                }else{
                    yPhi(6*j+ii) = 0;
                }
            }
        }

        DEInteg (VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
        //Extract state transition matrices
        for(int j=1;j<=6;j++){
            
            for(int i=1;i<=6;i++){
                Matrix aux = yPhi.slice(6*j+1,6*j+6);
                Phi(i,j) = aux(i);
            }
        }
        DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
        Y = Y_old;

         // Topocentric coordinates
        double theta = gmst(Mjd_UT1);                    // Earth rotation
        Matrix U = R_z(theta);
        Matrix r = Y.slice(1,3);
        Matrix s = LT*(U*r-Rs);                          // Topocentric position [m]

        // Time update
        TimeUpdate(P, Phi);
        
        // Azimuth and partials
        double Azim, Elev;
        Matrix dAds, dEds;
        AzElPa(Azim, Elev, dAds, dEds,s);     // Azimuth, Elevation
        Matrix dAdY(1,6);
        for(int i=1;i<=3;i++){
            dAdY(i) = (dAds*LT*U)(i);
        }
        for(int i=4;i<=6;i++){
            dAdY(i) = 0;
        }

        // Measurement update
        Matrix K;
        MeasUpdate (K, Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
        // Y.print();

        // Elevation and partials
        r = Y.slice(1,3);
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        AzElPa(Azim, Elev, dAds, dEds,s);     // Azimuth, Elevation
        
        Matrix dEdY(1,6);
        for(int i=1;i<=3;i++){
            dEdY(i) = (dEds*LT*U)(i);
        }
        for(int i=4;i<=6;i++){
            dEdY(i) = 0;
        }
        
        // Measurement update
        MeasUpdate ( K,Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );

        // Range and partials
        r = Y.slice(1,3);
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        double Dist = norm(s);
        Matrix dDds = ((1/Dist)*s).transponer();         // Range

        Matrix dDdY(1,6);
        for(int i=1;i<=3;i++){
            dDdY(i) = (dDds*LT*U)(i);
        }
        for(int i=4;i<=6;i++){
            dDdY(i) = 0;
        }
        
        // Measurement update
        MeasUpdate (K, Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );

    }
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,obs(46,1),'l');
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    double  Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);
    
    Matrix Y0 = Y;

    double Y_true_aux[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
    Matrix Y_true(6,1,Y_true_aux,6);
    cout << std::fixed << std::setprecision(1);
    printf("\nError of Position Estimation\n");
    cout << "dX      " << Y0(1)-Y_true(1)  <<  "[m]" << endl;
    cout << "dY      " << Y0(2)-Y_true(2)  <<  "[m]"  << endl;
    cout << "dZ      " << Y0(3)-Y_true(3)  <<  "[m]"  << endl;
    printf("\nError of Velocity Estimation\n");
    cout << "dVx      " << Y0(4)-Y_true(4)  <<  "[m/s]"  << endl;
    cout << "dVx      " << Y0(5)-Y_true(5)  <<  "[m/s]"  << endl;
    cout << "dVx      " << Y0(6)-Y_true(6)  <<  "[m/s]"  << endl;
}
