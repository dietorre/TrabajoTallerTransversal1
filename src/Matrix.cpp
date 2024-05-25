#include "../include/Matrix.h"

/**
 * @brief Default constructor.
 */
Matrix::Matrix() {}

/**
 * @brief Constructor with dimensions.
 * @param fil Number of rows.
 * @param col Number of columns.
 */
Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

/**
 * @brief Constructor with dimensions and initial values.
 * @param fil Number of rows.
 * @param col Number of columns.
 * @param v Array of initial values.
 * @param n Size of the array v.
 */
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

/**
 * @brief Copy constructor.
 * @param m Matrix to copy from.
 */
Matrix::Matrix(const Matrix& m)
{
    fil = m.fil;
    col = m.col;
    initMatrix();

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            matrix[i][j] = m.matrix[i][j];
        }
}

/**
 * @brief Constructor to concatenate two matrixes.
 * @param m1 First matrix.
 * @param m2 Second matrix.
 */
Matrix::Matrix(Matrix m1, Matrix m2){
    if(m1.col != m2.col){
        throw std::invalid_argument("Tamaños invalidos");
    }

    fil = m1.fil + m2.fil;
    col = m1.col;
    initMatrix();

    for(int i = 1; i <= m1.filas(); i++){
        for(int j = 1; j <= col; j++){
            matrix[i-1][j-1] = m1(i,j);
        }
    }
    
    for(int i = 1; i <= m2.filas(); i++){
        for(int j = 1; j <= col; j++){
            matrix[m1.filas() + i-1][j-1] = m2(i,j);
        }
    }
}

/**
 * @brief Constructor with an arithmetic sequence.
 * @param init Initial value.
 * @param salto Step size.
 * @param end End value.
 */
Matrix::Matrix(double init,double salto,double end)
{
    int nsaltos = floor((end-init)/salto);

    if(((end-init)/salto) - nsaltos != 0){
        nsaltos++;
    } 

    fil = nsaltos+1;
    col = 1;

    initMatrix();
    
    for(int i = 1; i <= nsaltos; i++){
        (*this)(i) = init + (i-1)*salto;
    }
    (*this)(nsaltos+1) = end;
}

/**
 * @brief Get the number of rows.
 * @return Number of rows.
 */
int Matrix::filas(){
    return fil;
}

/**
 * @brief Get the number of columns.
 * @return Number of columns.
 */
int Matrix::columnas(){
    return col;
}

/**
 * @brief Destructor.
 */
Matrix::~Matrix()
{   
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}

/**
 * @brief Initialize the matrix with zeros.
 */
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

/**
 * @brief Assignment operator.
 * @param matrix2 Matrix to assign from.
 * @return Reference to the assigned matrix.
 */
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    fil = matrix2.fil;
    col = matrix2.col;
    initMatrix();

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];
 
    return *this;
}

/**
 * @brief Addition operator.
 * @param matrix2 Matrix to add.
 * @return Result of the addition.
 */
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}

/**
 * @brief Addition operator with a scalar.
 * @param c Scalar to add.
 * @return Result of the addition.
 */
Matrix Matrix::operator+(const double c)
{
    Matrix result(fil,col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = c+matrix[i][j];

    return result;
}

/**
 * @brief Multiplication operator with a scalar.
 * @param c Scalar to multiply.
 * @param m2 Matrix to multiply.
 * @return Result of the multiplication.
 */
Matrix operator*(double c, Matrix m2)
{
    return m2*c;
}

/**
 * @brief Addition operator with a scalar.
 * @param c Scalar to add.
 * @param m2 Matrix to add.
 * @return Result of the addition.
 */
Matrix operator+(double c, Matrix m2)
{
    return m2+c;
}

/**
 * @brief Subtraction operator.
 * @param matrix2 Matrix to subtract.
 * @return Result of the subtraction.
 */
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}

/**
 * @brief Multiplication operator.
 * @param matrix2 Matrix to multiply.
 * @return Result of the multiplication.
 */
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, matrix2.col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}

/**
 * @brief Multiplication operator with a scalar.
 * @param c Scalar to multiply.
 * @return Result of the multiplication.
 */
Matrix Matrix::operator*(const double c)
{
    Matrix result(fil,col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = c*matrix[i][j];
 
    return result;
}

/**
 * @brief Get a reference to an element (2D).
 * @param i Row index.
 * @param j Column index.
 * @return Reference to the element.
 */
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}

/**
 * @brief Get a reference to an element (1D).
 * @param i Index.
 * @return Reference to the element.
 */
double& Matrix::operator()(const int i) const
{
    if(fil == 1){
        return matrix[0][i-1];
    }

    if(col == 1){
        return matrix[i-1][0];
    }

    throw std::invalid_argument("No es un vector");
}

/**
 * @brief Print the matrix.
 */
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief Get a specific row.
 * @param n Row index.
 * @return Row as a matrix.
 */
Matrix Matrix::getFila(int n){
    Matrix result(1, col);
    for(int i = 0; i < col; i++){
        result.matrix[0][i] = matrix[n-1][i]; 
    }
    return result;
}

/**
 * @brief Get a specific column.
 * @param n Column index.
 * @return Column as a matrix.
 */
Matrix Matrix::getColumna(int n){
    Matrix result(fil, 1);
    for(int i = 0; i < fil; i++){
        result.matrix[i][0] = matrix[i][n-1]; 
    }
    return result;
}

/** * @brief Calculate the Euclidean norm of the matrix (as a vector).
 * @return Euclidean norm.
 */
double Matrix::norm()
{
    if(col == 1){
        double sumaCuadrados = 0.0;
        for (size_t i = 0; i < fil; ++i) {
            sumaCuadrados += matrix[i][0] * matrix[i][0];
        }
        return sqrt(sumaCuadrados);
    }

    if(fil == 1){
        double sumaCuadrados = 0.0;
        for (size_t i = 0; i < col; ++i) {
            sumaCuadrados += matrix[0][i] * matrix[0][i];
        }
        return sqrt(sumaCuadrados);
    }
    throw std::invalid_argument("No es un vector");
}

/**
 * @brief Calculate the dot product of two vectors.
 * @param m First vector.
 * @param n Second vector.
 * @return Dot product.
 */
double dot(Matrix m, Matrix n){

    if(n.filas() != m.columnas()){
       throw std::invalid_argument("Tamaños inválidos");
    }

    if(m.filas() != 1 || n.columnas() != 1){
        throw std::invalid_argument("No son vectores");
    }

    double total = 0;

    for (int i = 1; i <= m.columnas(); i++) {
        total += m(i) * n(i);
    }

    return total;
}

/**
 * @brief Transpose the matrix.
 * @return Transposed matrix.
 */
Matrix Matrix::transponer(){

    Matrix res(col, fil);

    for (int i = 0; i < fil; i++)
    {
        for(int j = 0; j < col; j++){
            res.matrix[j][i] = matrix[i][j];
        }
    }

    return res;
}

/**
 * @brief Calculate the norm of a matrix (as a vector).
 * @param m Matrix to calculate the norm of.
 * @return Norm of the matrix.
 */
double norm(Matrix m){
    return m.norm();
}

/**
 * @brief Extract a slice of a matrix.
 * @param i Starting index.
 * @param j Ending index.
 * @return Sliced matrix.
 */
Matrix Matrix::slice(int i, int j){

    int tam;
    if(fil == 1){
        tam = col;
    }
    else if(col == 1){
        tam = fil;
    }
    else{
        throw std::invalid_argument("No es un vector");
    }

    if(i == -1){
        i = 1;
    }
    if(j == -1){
        j = tam;
    }

    Matrix res(j - i + 1, 1);

    for(int n = i; n <= j; n++){
        res(n - i + 1) = (*this)(n);
    } 

    return res;   
}

/**
 * @brief Calculate the inverse of the matrix.
 * @return Inverse of the matrix.
 */
Matrix Matrix::inverse() {
    if (fil != col) {
        throw std::runtime_error("La matriz no es cuadrada, no se puede calcular la inversa.");
    }

    int n = fil;

    // Create an extended matrix with the original matrix and the identity matrix
    Matrix extended(n, 2 * n);
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            extended(i, j) = matrix[i-1][j-1];
        }
        extended(i, i + n) = 1; // Fill the right part with the identity matrix
    }

    // Apply Gaussian elimination to obtain the reduced row-echelon form
    for (int i = 1; i <= n; ++i) {
        // Make the diagonal element of the i-th submatrix equal to 1
        double pivot = extended(i, i);
        if (pivot == 0) {
            throw std::runtime_error("La matriz es singular, no se puede calcular la inversa.");
        }
        for (int j = 1; j <= 2 * n; ++j) {
            extended(i, j) /= pivot;
        }

        // Make zeros in the i-th column
        for (int k = 1; k <= n; ++k) {
            if (k != i) {
                double factor = extended(k, i);
                for (int j = 1; j <= 2 * n; ++j) {
                    extended(k, j) -= factor * extended(i, j);
                }
            }
        }
    }

    // Extract the inverse matrix from the right part of the extended matrix
    Matrix result(n, n);
    for (int i = 1; i <= n; ++i) {
        for (int j = 1;j <= n; ++j) {
            result(i, j) = extended(i, j + n);
        }
    }

    return result;
}

