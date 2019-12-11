/**
 * @file vecMath3D
 * @brief Operator overloads and functions used for vector analysis in three dimensions.
 *
 * Contains addition and subtracting of vectors with other vectors and scalars,
 * matrix and vector multiplication, calculating the absolute value of vectors and
 * multiplication of matrices with vectors.
 */
#include "vecMath.hpp"

typedef std::array<double, 3> Vector;
typedef std::array<Vector, 3> Matrix;

Vector operator+ (Vector vecA, Vector vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1], vecA[2] + vecB[2] };
}

array<int, 3> operator+ (array<int, 3> vecA, array<int, 3> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1], vecA[2] + vecB[2] };
}

Vector operator- (Vector vecA, Vector vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1], vecA[2] - vecB[2] };
}

Vector operator- (Vector vecA, array<int, 3> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1], vecA[2] - vecB[2] };
}

Vector operator* (double a, const Vector& vec) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

Vector operator* (const Vector& vec, double a) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

Vector operator* (const array<int, 3>& vec, double a) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

array<Vector, 3> operator* (double a, const array<array<double,3>, 3>& matrix) {
    array<Vector, 3> tempReturn;
	for (int row = 0; row < 3; ++row)
		for (int col = 0; col < 3; ++col)
			tempReturn[row][col] = a*matrix[row][col];

	return tempReturn;
}

double operator* (const Vector& vecA, const Vector& vecB) {
  return vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2];
}

Vector operator* (const array<array<double,3>, 3>& matrix, const Vector vec) {
    Vector tempReturn;
    for (int row = 0; row < 3; row++)
		tempReturn[row] = matrix[row]*vec;

    return tempReturn;
}

Matrix operator* (const array<array<double,3>, 3>& matrixA, const array<array<double,3>, 3>& matrixB) {
    Matrix tempReturn;
    Matrix matrixB_T = transpose(matrixB);

    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
            tempReturn[row][col] = matrixA[row]*matrixB_T[col];

    return tempReturn;
}

Vector operator/ (const Vector& vec, double a) {
  return { vec[0]/a, vec[1]/a, vec[2]/a };
}

double abs (const Vector vec) {
  return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}

array<Vector, 3> transpose(array<Vector, 3> matrix) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < i; j++)
			std::swap(matrix[i][j], matrix[j][i]);

	return matrix;
}
