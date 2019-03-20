#include "vecMath3D.hpp"

std::array<double, 3> operator+ (std::array<double, 3> vecA, std::array<double, 3> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1], vecA[2] + vecB[2] };
}

std::array<int, 3> operator+ (std::array<int, 3> vecA, std::array<int, 3> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1], vecA[2] + vecB[2] };
}

std::array<double, 3> operator- (std::array<double, 3> vecA, std::array<double, 3> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1], vecA[2] - vecB[2] };
}

std::array<double, 3> operator- (std::array<double, 3> vecA, std::array<int, 3> vecB){
  return { vecA[0] - vecB[0], vecA[1] - vecB[1], vecA[2] - vecB[2] };
}

std::array<double, 3> operator* (double a, const std::array<double, 3>& vec) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

std::array<double, 3> operator* (const std::array<double, 3>& vec, double a) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

std::array<double, 3> operator* (const std::array<int, 3>& vec, double a) {
  return { a*vec[0], a*vec[1], a*vec[2] };
}

double operator* (const std::array<double, 3>& vecA, const std::array<double, 3>& vecB) {
  return { vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2] };
}

std::array<double, 3> operator/ (const std::array<double, 3>& vec, double a) {
  return { vec[0]/a, vec[1]/a, vec[2]/a };
}

double abs (const std::array<double, 3> vec) {
  return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}
