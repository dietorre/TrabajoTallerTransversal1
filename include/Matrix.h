#ifndef _MATRIX_
#define _MATRIX_

#include <cmath>
#include <iostream>
#include <iomanip>

class Matrix
{
    public:
        Matrix();
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        Matrix(const Matrix& m);
        Matrix(Matrix m1, Matrix m2);
        Matrix(double init,double salto,double end);
        ~Matrix();
        int filas(); 
        int columnas(); 
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator+(const double c);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        Matrix  operator*(const double c);
        friend Matrix operator*(double c, Matrix m2);
        friend Matrix operator+(double c, Matrix m2);
        double& operator()(const int i, const int j) const;
        double& operator()(const int i) const;
        

        Matrix getFila(int n);
        Matrix getColumna(int n);
        double norm();
        Matrix transponer();
        Matrix slice(int i, int j);
        Matrix inverse();
 
        void print();
 
    private:
        void initMatrix();
        int fil;
        int col;
        double **matrix;
};

double dot(Matrix m, Matrix n);
double norm(Matrix m);

#endif
