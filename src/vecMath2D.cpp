#include "vecMath2D.h"

std::array<double, 2> operator+ (std::array<double, 2> vecA, std::array<double, 2> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1]};
}

std::array<double, 2> operator- (std::array<double, 2> vecA, std::array<double, 2> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1]};
}

std::array<double, 2> operator- (std::array<double, 2> vecA, std::array<int, 2> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1]};
}

std::array<double, 2> operator*(double a, const std::array<double, 2>& vec) {
  return { a*vec[0], a*vec[1]};
}

std::array<double, 2> operator*(const std::array<double, 2>& vec, double a) {
  return { a*vec[0], a*vec[1]};
}

std::array<double, 2> operator*(const std::array<int, 2>& vec, double a) {
  return { a*vec[0], a*vec[1]};
}

double operator* (const std::array<double, 2>& vecA, const std::array<double, 2>& vecB) {
  return { vecA[0]*vecB[0] + vecA[1]*vecB[1]};
}

double abs (const std::array<double, 2> vec) {
  return sqrt(pow(vec[0], 2) + pow(vec[1], 2));
}
