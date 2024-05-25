#include "../include/initGlobales.h"

using namespace std;

void eop19620101()
{
    extern Matrix eopdata;
    FILE *fp = fopen("./data/eop19620101.txt","r");
    if(fp == NULL){
        printf("Error");
        exit(EXIT_FAILURE);
    }

    for (int i=1;i<=21413;i++){
        fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(eopdata(1,i)),
               &(eopdata(2,i)),&(eopdata(3,i)),&(eopdata(4,i)),
               &(eopdata(5,i)),&(eopdata(6,i)),&(eopdata(7,i)),
               &(eopdata(8,i)),&(eopdata(9,i)),&(eopdata(10,i)),
               &(eopdata(11,i)),&(eopdata(12,i)),&(eopdata(13,i)));
    }

    fclose(fp);
}

void GGM03S()
{
    extern Matrix Cnm;
    extern Matrix Snm;
    FILE *fp = fopen("./data/GGM03S.txt","r");
    if(fp == NULL){
        printf("Error");
        exit(EXIT_FAILURE);
    }

    int aux1,aux2;
    double aux3,aux4;

    for(int n=0;n<=180;n++){
        for(int m=0;m<=n;m++){
            fscanf(fp,"%d %d %lf %lf %lf %lf",&aux1,&aux2,&(Cnm(n+1,m+1)),&(Snm(n+1,m+1)),&aux3,&aux4);
        }
    }

    fclose(fp);
}

void DE430Coeff()
{
    extern Matrix PC;

    ifstream file;

    file.open("./data/M_tab.txt");

    for(int i=1;i<=2285;i++){
        for(int j=1;j<=1020;j++){
            file >> PC(i,j);
        }
    }
    file.close();
}

void GEOS3(){

    extern int nobs;
    extern Matrix obs;

    FILE *fp = fopen("./data/GEOS3.txt","r");

    if(fp == NULL){
        printf("Error");
        exit(EXIT_FAILURE);
    }

    char *tline = new char[200];
    fgets(tline,200,fp);

    char* out;
    int Y,M,D,h,m;
    double s,az,el,Dist;

    for(int i=1;i<=nobs;i++){
        char* nueva = new char[10];

        strncpy(nueva,tline,4);
        nueva[4] = '\0';
        Y = strtol(nueva,&out,10);

        strncpy(nueva,tline+5,2);
        nueva[2] = '\0';
        M = strtol(nueva,&out,10);

        strncpy(nueva,tline+8,2);
        nueva[2] = '\0';
        D = strtol(nueva,&out,10);

        strncpy(nueva,tline+12,2);
        nueva[2] = '\0';
        h = strtol(nueva,&out,10);

        strncpy(nueva,tline+15,2);
        nueva[2] = '\0';
        m = strtol(nueva,&out,10);

        strncpy(nueva,tline+18,5);
        nueva[5] = '\0';
        s = atof(nueva);

        strncpy(nueva,tline+25,8);
        nueva[8] = '\0';
        az = atof(nueva);

        strncpy(nueva,tline+35,8);
        nueva[8] = '\0';
        el = atof(nueva);

        strncpy(nueva,tline+44,12);
        nueva[12] = '\0';
        Dist = atof(nueva);

        obs(i,1) = mjday(Y,M,D,h,m,s);
        obs(i,2) = Const::Rad*az;
        obs(i,3) = Const::Rad*el;
        obs(i,4) = 1e3*Dist;

        fgets(tline,100,fp);
    }

    fclose(fp);
}
