#include "include/Matrix.h"
#include "include/R_z.h"
#include "include/R_x.h"
#include "include/R_y.h"
#include "include/Mjday.h"
#include "include/Mjday_TDB.h"
#include "include/Frac.h"
#include "include/sign_.h"
#include "include/SAT_Const.h"
#include "include/IERS.h"
#include "include/AccelPointMass.h"
#include "include/EccAnom.h"
#include "include/Position.h"
#include "include/AzElPa.h"
#include "include/Cheb3D.h"
#include "include/Geodetic.h"
#include "include/MeanObliquity.h"
#include "include/Legendre.h"
#include "include/NutAngles.h"
#include "include/unit.h"
#include "include/globales.h"
#include "include/PrecMatrix.h"
#include "include/NutMatrix.h"
#include "include/AccelHarmonic.h"
#include "include/gmst.h"
#include "include/EqnEquinox.h"
#include "include/gast.h"
#include "include/GHAMatrix.h"
#include "include/JPL_Eph_DE430.h"
#include "include/Accel.h"
#include "include/LTC.h"
#include "include/G_AccelHarmonic.h"
#include "include/VarEqn.h"
#include "include/MeasUpdate.h"
#include "include/initGlobales.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)
#define TOL_ 10e-14
using namespace std;

int test_Matrix_constructor() {
    // Crear una matriz 3x3 con valores específicos
    double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix mat(3, 3, values, 9);

    // Valores esperados
    double expected_values[3][3] = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };

    // Verificar que cada elemento de la matriz coincide con el esperado
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            _assert(fabs(mat(i+1, j+1) - expected_values[i][j]) < TOL_);
        }
    }

    return 0;
}

int test_Matrix_multiplication() {
    // Crear dos matrices 2x2
    double values1[] = {1.0, 2.0, 3.0, 4.0};
    double values2[] = {5.0, 6.0, 7.0, 8.0};
    Matrix mat1(2, 2, values1, 4);
    Matrix mat2(2, 2, values2, 4);

    // Realizar la multiplicación
    Matrix result = mat1 * mat2;

    // Valores esperados
    double expected_values[2][2] = {
        {19.0, 22.0},
        {43.0, 50.0}
    };

    // Verificar que cada elemento de la matriz resultante coincide con el esperado
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            _assert(fabs(result(i+1, j+1) - expected_values[i][j]) < TOL_);
        }
    }

    return 0;
}

int test_Matrix_transpose() {
    // Crear una matriz 2x3
    double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    Matrix mat(2, 3, values, 6);

    // Realizar la transposición
    Matrix transposed = mat.transponer();

    // Valores esperados
    double expected_values[3][2] = {
        {1.0, 4.0},
        {2.0, 5.0},
        {3.0, 6.0}
    };

    // Verificar que cada elemento de la matriz transpuesta coincide con el esperado
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            _assert(fabs(transposed(i+1, j+1) - expected_values[i][j]) < TOL_);
        }
    }

    return 0;
}

int test_Matrix_norm() {
    // Crear una matriz 3x1 con valores específicos
    double values[] = {3.0, 4.0, 5.0};
    Matrix mat(3, 1, values, 3);

    // Calcular la norma
    double result = mat.norm();

    // Valor esperado (norma Euclidiana)
    double expected = sqrt(3.0*3.0 + 4.0*4.0 + 5.0*5.0);

    // Verificar que el resultado coincide con el valor esperado
    _assert(fabs(result - expected) < TOL_);

    return 0;
}



int test_R_z_1(){
    double angle = 2.0;
    Matrix sol(3,3);
    sol = R_z(angle);
    //sol.print();

    _assert((fabs(sol(1,1) + 0.416146836547142) < TOL_) && (fabs(sol(1,2) - 0.909297426825682) < TOL_) && (fabs(sol(1,3)) < TOL_));
    _assert((fabs(sol(2,1) + 0.909297426825682) < TOL_) && (fabs(sol(2,2) + 0.416146836547142) < TOL_) && (fabs(sol(2,3)) < TOL_));
    _assert((fabs(sol(3,1)) < TOL_) && (fabs(sol(3,2)) < TOL_) && (fabs(sol(3,3) - 1.0) < TOL_));

    return 0;
}

int test_R_x_1(){
    double angle = 2.0;
    Matrix sol(3,3);
    sol = R_x(angle);

    _assert((fabs(sol(1,1) - 1) < TOL_) && (fabs(sol(1,2)) < TOL_) && (fabs(sol(1,3)) < TOL_));
    _assert((fabs(sol(2,1)) < TOL_) && (fabs(sol(2,2) + 0.416146836547142) < TOL_) && (fabs(sol(2,3) - 0.909297426825682)  < TOL_));
    _assert((fabs(sol(3,1)) < TOL_) && (fabs(sol(3,2) + 0.909297426825682)  < TOL_) && (fabs(sol(3,3) + 0.416146836547142)  < TOL_));

    return 0;
}

int test_R_y_1(){
    double angle = 2.0;
    Matrix sol(3,3);
    sol = R_y(angle);

    _assert((fabs(sol(1,1) + 0.416146836547142) < TOL_) && (fabs(sol(1,2)) < TOL_) && (fabs(sol(1,3) + 0.909297426825682)  < TOL_));
    _assert((fabs(sol(2,1)) < TOL_) && (fabs(sol(2,2) - 1) < TOL_) && (fabs(sol(2,3))  < TOL_));
    _assert((fabs(sol(3,1) - 0.909297426825682) < TOL_) && (fabs(sol(3,2))  < TOL_) && (fabs(sol(3,3) + 0.416146836547142)  < TOL_));

    return 0;
}

int test_frac(){
    double x = 1.23456789;
    _assert(fabs(Frac(x) - 0.23456789) < TOL_);

    return 0;
}

int test_sign_1(){
    _assert(fabs(sign_(-1,-1) + 1) < TOL_);
    return 0;
}

int test_sign_2(){
    _assert(fabs(sign_(-1,1) - 1) < TOL_);
    return 0;
}

int test_sign_3(){
    _assert(fabs(sign_(1,-1) + 1) < TOL_);
    return 0;
}

int test_sign_4(){
    _assert(fabs(sign_(1,1) - 1) < TOL_);
    return 0;
}

int test_mjday() {
    _assert(fabs(mjday(2024, 3, 23) - 60392) <= TOL_);
    return 0;
}

int test_Mjday_TDB() {
//    cout << setprecision(14);
//    cout << Mjday_TDB(1);
    _assert(fabs(Mjday_TDB(1) - 0.999999986586878) <= TOL_);
    return 0;
}

int test_mjday_2() {
    _assert(fabs(mjday(2024, 02, 29, 15, 12, 12) - 6.036963347222237e+04) <= TOL_);
    return 0;
}

int test_mjday_3() {
    _assert(fabs(mjday(1995, 1, 29, 02, 38, 0) - 4.974610972222220e+04) <= TOL_);
    return 0;
}

int test_IERS(){

    extern Matrix eopdata;
    double x_pole = 0;
    double y_pole= 0;
    double UT1_UTC= 0;
    double LOD= 0;
    double dpsi= 0;
    double deps= 0;
    double dx_pole= 0;
    double dy_pole= 0;
    double TAI_UTC= 0;

    double Mjd_UTC = 49746.1163541665;

    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,Mjd_UTC,'l');

    _assert(fabs(x_pole + 5.5937872420407e-07) < TOL_);
    _assert(fabs(y_pole - 2.33559834147197e-06) < TOL_);
    _assert(fabs(UT1_UTC - 0.325747632958709) < TOL_);
    _assert(fabs(LOD - 0.00272698971874332) < TOL_);
    _assert(fabs(dpsi + 1.16882953161744e-07) < TOL_);
    _assert(fabs(deps + 2.4783506198648e-08) < TOL_);
    _assert(fabs(dx_pole + 8.43027359626024e-10) < TOL_);
    _assert(fabs(dy_pole + 1.56811369105037e-09) < TOL_);
    _assert(fabs(TAI_UTC -29) < TOL_);
    return 0;
}

int test_norm(){
    Matrix m(3,1);
    m(1,1) = 1;
    m(2,1) = 2;
    m(3,1) = 3;
    _assert(fabs(m.norm() -  3.74165738677394) < TOL_);
    return 0;
}

int test_accel_point_mass(){
    Matrix r(1,3);
    r(1,1) = 7; r(1,2) = 8; r(1,3) = 9;
    
    Matrix s(1,3);
    s(1,1) = 3; s(1,2) = 2; s(1,3) = 1;

    Matrix v = AccelPointMass(r,s,10.0);

    _assert(
        fabs(v(1,1) + 0.604719098857642) < TOL_ &&
        fabs(v(1,2) + 0.429826430585706) < TOL_ &&
        fabs(v(1,3) + 0.254933762313769) < TOL_
    );
    
    return 0;
}

int test_eccanom(){
    _assert(fabs(EccAnom(15.5,17.0) - 3.13003889593541) < TOL_);
    return 0;
}

int test_Position(){
    Matrix p = Position(100,2000,50);
    _assert(fabs(p(1) + 2026915.55906022922135889530181884765625) < TOL_ &&
            fabs(p(2) - 1190233.021128253079950809478759765625) < TOL_ &&
            fabs(p(3) - 5909388.454115637578070163726806640625) < TOL_);
    return 0;
}

int test_AzElPa(){
    double Az;
    double El;
    Matrix dAds(3,1);
    Matrix dEds(3,1);
    Matrix s(3,1);
    s(1)=1;s(2)=2;s(3)=3;
    AzElPa(Az,El,dAds,dEds,s);
    _assert(fabs(Az - 0.463647609000806) < TOL_);
    _assert(fabs(El - 0.930274014115472) < TOL_);
    _assert(fabs(dAds(1) -  0.4) < TOL_ && fabs(dAds(2) + 0.2) < TOL_ && fabs(dAds(3)) < TOL_);
    _assert(fabs(dEds(1) + 0.095831484749991) < TOL_ && fabs(dEds(2) + 0.191662969499982) < TOL_ && fabs(dEds(3) - 0.159719141249985) < TOL_);
    return 0;
}

int test_Cheb(){
    double t = 0.5;
    double N = 3;
    double Ta = 0;
    double Tb = 1;
    Matrix Cx(3,1);
    Cx(1)=10;Cx(2)=11;Cx(3)=12;

    Matrix Cy(3,1);
    Cy(1)=-1;Cy(2)=-2;Cy(3)=-3;

    Matrix Cz(3,1);
    Cz(1)=0;Cz(2)=9;Cz(3)=-9;

    Matrix res = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);
    _assert(fabs(res(1) + 2) < TOL_);
    _assert(fabs(res(2) - 2) < TOL_);
    _assert(fabs(res(3) - 9) < TOL_);
    return 0;
}

int test_geodetic(){
    double lon,lat,h;
    Matrix r(3,1);
    r(1)=7;r(2)=8;r(3)=9;
    Geodetic(lon,lat,h,r);
    _assert(fabs(lon - 0.851966327173272) < TOL_);
    _assert(fabs(lat - 1.57054825048658) < TOL_);
    _assert(fabs(h + 6356742.615273603238165378570556640625) < TOL_);
    return 0;
}

int test_MeanObliquity(){
    _assert(fabs(MeanObliquity(100) - 0.409412448949153) < TOL_);
    return 0;
}

int test_Legendre(){
    Matrix pnm(3,3);
    Matrix dpnm(3,3);
    Legendre(pnm,dpnm,2,2,10);

    _assert(fabs(pnm(1,1) - 1) < TOL_ && fabs(pnm(1,2)) < TOL_ && fabs(pnm(1,3)) < TOL_);
    _assert(fabs(pnm(2,1) + 0.942272204450451) < TOL_ && fabs(pnm(2,2) + 1.45331451954492) < TOL_ && fabs(pnm(2,3)) < TOL_);
    _assert(fabs(pnm(3,1) + 0.125357428584815) < TOL_ && fabs(pnm(3,2) - 1.76791087603363) < TOL_ && fabs(pnm(3,3) - 1.36336959387417) < TOL_);

    _assert(fabs(dpnm(1,1)) < TOL_ && fabs(dpnm(1,2)) < TOL_ && fabs(dpnm(1,3)) < TOL_);
    _assert(fabs(dpnm(2,1) + 1.45331451954492) < TOL_ && fabs(dpnm(2,2) - 0.942272204450451) < TOL_ && fabs(dpnm(2,3)) < TOL_);
    _assert(fabs(dpnm(3,1) - 3.06211146054385) < TOL_ && fabs(dpnm(3,2) - 1.58049502928925) < TOL_ && fabs(dpnm(3,3) + 1.76791087603363) < TOL_);

    // pnm.print();
    // dpnm.print();
    return 0;
}

int test_NutAngles(){
    double dpsi,deps;
    NutAngles(dpsi,deps,15);
    _assert(fabs( dpsi - 3.07180427954353e-05 ) < TOL_);
    _assert(fabs( deps - 3.78674050874946e-05 ) < TOL_);
    return 0;
}

int test_unit(){
    Matrix vec(3,1);
    vec(1,1)=4;
    vec(2,1)=5;
    vec(3,1)=6;
    Matrix res = unit(vec);
    _assert(fabs(res(1) - 0.455842305838552) < TOL_);
    _assert(fabs(res(2) - 0.56980288229819) < TOL_);
    _assert(fabs(res(3) - 0.683763458757828) < TOL_);
    return 0;
}

int test_PrecMatrix(){
    Matrix res(3,3);
    res = PrecMatrix(10,5);
    _assert(fabs(res(1,1) - 0.999999999994437) < TOL_ && fabs(res(1,2) - 3.05853689370122e-06) < TOL_ && fabs(res(1,3) - 1.33100733379074e-06) < TOL_);
    _assert(fabs(res(2,1) + 3.05853689370122e-06) < TOL_ && fabs(res(2,2) - 0.999999999995323) < TOL_ && fabs(res(2,3) + 2.03546747019913e-12) < TOL_);
    _assert(fabs(res(3,1) + 1.33100733379074e-06) < TOL_ && fabs(res(3,2) + 2.03546756599807e-12) < TOL_ && fabs(res(3,3) - 0.999999999999114) < TOL_);
    return 0;
}

int test_NutMatrix(){
    Matrix res(3,3);
    res = NutMatrix(10);
    _assert(fabs(res(1,1) - 0.999999999492159) < TOL_ && fabs(res(1,2) + 2.92358797494727e-05) < TOL_ && fabs(res(1,3) + 1.26864277798066e-05) < TOL_);
    _assert(fabs(res(2,1) - 2.9235393806329e-05) < TOL_ && fabs(res(2,2) - 0.999999998839099) < TOL_ && fabs(res(2,3) + 3.83026695232602e-05) < TOL_);
    _assert(fabs(res(3,1) - 1.26875475773192e-05) < TOL_ && fabs(res(3,2) - 3.83022986110704e-05) < TOL_ && fabs(res(3,3) - 0.99999999918598) < TOL_);
    return 0;
}

int test_AccelHarmonic(){
    Matrix a(3,1),r(3,1),E(3,3);
    
    r(1) = 5720694.2260585; r(2) = 2687728.41425142; r(3) = 3483000.08675422;

    E(1,1) = -0.976558757940107; E(1,2) = 0.215250556888025; E(1,3) = -0.000435947096290288;
    E(2,1) = -0.2152505181354; E(2,2) = -0.976558854525697; E(2,3) = -0.000134498699131133;
    E(3,1) = -0.000454678916875739; E(3,2) = -3.7508044211961e-05; E(3,3) = 0.999999895930109;

    a = AccelHarmonic(r,E,20,20);
    
    _assert(fabs(a(1) + 6.06544113186609) < TOL_);
    _assert(fabs(a(2) + 2.8497772091355) < TOL_);
    _assert(fabs(a(3) + 3.70232504408661) < TOL_);
    
    return 0;
}

int test_gmst(){
    _assert(fabs(gmst(15) - 1.23125001740142) < TOL_);
    return 0;
}

int test_EqnEquinox(){
    _assert(fabs(EqnEquinox(-1) - 2.59193273800512e-05) < TOL_);
    return 0;
}

int test_gast(){
    _assert(fabs(gast(-1) - 0.956031276214774) < TOL_);
    return 0;
}

int test_GHAMatrix(){
    
    Matrix res = GHAMatrix(49746.11006769535015337169170379638671875);
    _assert(fabs(res(1,1) + 0.976333634522934) < TOL_);
    _assert(fabs(res(1,1) + 0.976333634522934) < TOL_ && fabs(res(1,2) - 0.216269817818478) < TOL_ && fabs(res(1,3)) < TOL_);
    _assert(fabs(res(2,1) + 0.216269817818478) < TOL_ && fabs(res(2,2) + 0.976333634522934) < TOL_ && fabs(res(2,3)) < TOL_);
    _assert(fabs(res(3,1)) < TOL_ && fabs(res(3,2)) < TOL_ && fabs(res(3,3) - 1) < TOL_);
    return 0;
}

int test_JPL_Eph(){
    Matrix r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun;
    JPL_Eph_DE430(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun,49746.1107720813088235445320606231689453125);
    cout << setprecision(100);
    // _assert(fabs(r_Mercury(1) - 83779075648.5296783447265625) < TOL_);
    //   && fabs(r_Mercury(2) + 65292051794.8271636962890625) < TOL_ && fabs(r_Mercury(2) + 23392255339.134342193603515625) < TOL_);
    // r_Mercury.print();
    // r_Venus.print();
    // r_Earth.print();
    // r_Mars.print();
    // r_Jupiter.print();
    // r_Saturn.print();
    // r_Uranus.print();
    // r_Neptune.print();
    // r_Pluto.print();
    // r_Moon.print();
    // r_Sun.print();
    return 0;
}

int test_Accel(){
    Matrix y(6,1);
    y(1) = 6221397.62857869;
    y(2) = 2867713.77965738;
    y(3) = 3006155.98509949;
    y(4) = 4645.04725161806;
    y(5) = -2752.21591588204;
    y(6) =  -7507.99940987031;

    Matrix res = Accel(0,y);
    _assert(fabs(res(1) - 4645.04725161806) < TOL_);
    _assert(fabs(res(2) + 2752.21591588204) < TOL_);
    _assert(fabs(res(3) + 7507.99940987031) < TOL_);
    _assert(fabs(res(4) + 5.92414951314006) < TOL_);
    _assert(fabs(res(5) + 2.73076669788113) < TOL_);
    _assert(fabs(res(6) + 2.86933570556259) < TOL_);
    
    return 0;
}

int test_LTC(){
    Matrix res = LTC(-2.76234307910694,0.376551295459273);
    _assert(fabs(res(1,1) - 0.370223471399199) < TOL_ && fabs(res(1,2) + 0.928942722252092) < TOL_ && fabs(res(1,3)) < TOL_);
    _assert(fabs(res(2,1) - 0.341586711932422) < TOL_ && fabs(res(2,2) - 0.136136938528208) < TOL_ && fabs(res(2,3) - 0.929938305587722) < TOL_);
    _assert(fabs(res(3,1) + 0.863859421119156) < TOL_ && fabs(res(3,2) + 0.344284987681776) < TOL_ && fabs(res(3,3) - 0.367715580035218) < TOL_);
    return 0;
}

int test_G_AccelHarmonic(){
    Matrix r(3,1),U(3,3);
    r(1) = 7101800.90695316;
    r(2) = 1293997.58115302;
    r(3) = 10114.0149489548;

    U(1,1)=-0.984320311904791; U(1,2)=0.17638970840918; U(1,3)=-0.000440838949610109;
    U(2,1)=-0.176389673507182; U(2,2)=-0.984320409906027; U(2,3)=-0.000117142904888635;
    U(3,1)=-0.000454589578418276; U(3,2)=-3.75467022865179e-05; U(3,3)=0.999999895969275;

    Matrix res = G_AccelHarmonic(r,U,20,20);
    
    _assert(fabs(res(1,1) - 2.02233500345983e-06) < TOL_ && fabs(res(1,2) - 5.61803299881092e-07) < TOL_ && fabs(res(1,3) - 4.3985650677314e-09) < TOL_);
    _assert(fabs(res(2,1) - 5.61803300325181e-07) < TOL_ && fabs(res(2,2) + 9.58631633629636e-07) < TOL_ && fabs(res(2,3) - 8.05635336220689e-10) < TOL_);
    _assert(fabs(res(3,1) - 4.39855909334375e-09) < TOL_ && fabs(res(3,2) - 8.05634045586423e-10) < TOL_ && fabs(res(3,3) + 1.06370336962376e-06) < TOL_);

    return 0;
}


int all_tests() {
    _verify(test_mjday);
    _verify(test_mjday_2);
    _verify(test_mjday_3);
    _verify(test_R_z_1);
    _verify(test_R_x_1);
    _verify(test_R_y_1);
    _verify(test_frac);
    _verify(test_sign_1);
    _verify(test_sign_2);
    _verify(test_sign_3);
    _verify(test_sign_4);
    _verify(test_Mjday_TDB);
    _verify(test_IERS);
    _verify(test_norm);
    _verify(test_accel_point_mass);
    _verify(test_eccanom);
    _verify(test_Position);
    _verify(test_AzElPa);
    _verify(test_Cheb);
    _verify(test_geodetic);
    _verify(test_MeanObliquity);
    _verify(test_Legendre);
    _verify(test_NutAngles);
    _verify(test_unit);
    _verify(test_PrecMatrix);
    _verify(test_NutMatrix);
    _verify(test_AccelHarmonic);
    _verify(test_gmst);
    _verify(test_EqnEquinox);
    _verify(test_gast);
    _verify(test_GHAMatrix);
    _verify(test_JPL_Eph);
    _verify(test_Accel);
    _verify(test_LTC);
    _verify(test_G_AccelHarmonic);
    _verify(test_Matrix_constructor);
    _verify(test_Matrix_multiplication);
    _verify(test_Matrix_transpose);
    _verify(test_Matrix_norm);
    return 0;
}

int main() {

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

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    /*Global::eop19620101(21413);
    Global::GGM03S();
    Global::GEOS3();*/
    // obs.print();
    //(*Global::obs).print();
    //Global::DE430Coeff();
/*
    FILE *fp = fopen("./data/GEOS3.txt","r");
    if(fp == NULL){
        printf("Error");
        exit(EXIT_FAILURE);
    }
    char *tline = new char[100];

    fgets(tline,100,fp);
    fclose(fp);

    cout << tline;
//    int tam = 2;
    char* nueva = new char[10];
    strncpy(nueva,tline,4);
    nueva[4] = '\0';
    char* out;
    cout << strtol(nueva,&out,10) << endl;*/

    // Matrix m(3,3);
    // m(1,1) = 1; m(1,2) = 1; m(1,3) = 0;
    // m(2,1) = 1; m(2,2) = 0; m(2,3) = 1;
    // m(3,1) = 0; m(3,2) = 1; m(3,3) = 0;

    // m.inverse().print();
    
    // cout << m.getFila(1).norm() << endl;

    // Matrix temp(171,10,201);
    // temp.print();
    
    // Matrix r(5,1);
    // r(1) = 1; r(2) = 2; r(3) = 3; r(4) = 4; r(5) = 5; 

    // Matrix s(2,1);
    // s(1)=15;s(2)=10;

    // Matrix aux(s,r);

    // aux.print();
    // r.slice(5,-1).print();

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}