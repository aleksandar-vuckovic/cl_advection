#ifndef VEC_MATH_3D
#define VEC_MATH_3D

#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

using std::array;

typedef array<double, 3> Vector;
typedef array<Vector, 3> Matrix;

Vector operator+ (Vector vecA, Vector vecB);

array<int, 3> operator+ (array<int, 3> vecA, array<int, 3> vecB);

Vector operator- (Vector vecA, Vector vecB);

Vector operator- (Vector vecA, array<int, 3> vecB);

Vector operator*(double a, const Vector& vec);

Vector operator*(const Vector& vec, double a);

Vector operator*(const array<int, 3>& vec, double a);

Matrix operator* (double a, const array<array<double,3>, 3>& matrix);

double operator* (const Vector& vecA, const Vector& vecB);

Vector operator* (const array<array<double,3>, 3>& matrix, const Vector vec);

Matrix operator* (const Matrix& matrixA, const Matrix& matrixB);

Vector operator/ (const Vector& vec, double a);

double abs (const Vector vec);

Matrix transpose(const Matrix matrix);

#endif
