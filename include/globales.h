#ifndef TRABAJOTT1_GLOBALES_H
#define TRABAJOTT1_GLOBALES_H

#include "Matrix.h"
#include "SAT_Const.h"
#include "Mjday.h"
#include "parametros.h"

#include <cstdlib>
#include <fstream>
#include <cstring>

Matrix eopdata(13,21413);
Matrix Cnm(181,181);
Matrix Snm(181,181);
Matrix PC(2285,1020);
int nobs = 46;
Matrix obs(nobs,4);

Parametros AuxParam;

#endif //TRABAJOTT1_GLOBALES_H
