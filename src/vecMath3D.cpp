#include "vecMath3D.hpp"

using std::array;

array<double, 3> operator+ (array<double, 3> vecA, array<double, 3> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1], vecA[2] + vecB[2] };
}

array<int, 3> operator+ (array<int, 3> vecA, array<int, 3> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1], vecA[2] + vecB[2] };
}

array<double, 3> operator- (array<double, 3> vecA, array<double, 3> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1], vecA[2] - vecB[2] };
}

array<double, 3> operator- (array<double, 3> vecA, array<int, 3> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1], vecA[2] - vecB[2] };
}

array<double, 3> operator* (double a, const array<double, 3>& vec) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

array<double, 3> operator* (const array<double, 3>& vec, double a) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

array<double, 3> operator* (const array<int, 3>& vec, double a) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

array<array<double, 3>, 3> operator* (double a, const array<array<double,3>, 3>& matrix) {
	array<array<double, 3>, 3> tempReturn;
	for (int row = 0; row < 3; ++row)
		for (int col = 0; col < 3; ++col)
			tempReturn[row][col] = a*matrix[row][col];

	return tempReturn;
}

double operator* (const array<double, 3>& vecA, const array<double, 3>& vecB) {
  return { vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2] };
}

array<double, 3> operator* (const array<array<double,3>, 3>& matrix, const array<double, 3> vec) {
    array<double, 3> tempReturn;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
            tempReturn[row] = matrix[row]*vec;

    return tempReturn;
}

array<double, 3> operator/ (const array<double, 3>& vec, double a) {
  return { vec[0]/a, vec[1]/a, vec[2]/a };
}

double abs (const array<double, 3> vec) {
  return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}
