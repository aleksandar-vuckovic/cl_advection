#ifndef VEC_MATH_3D
#define VEC_MATH_3D

#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

std::array<double, 3> operator+ (std::array<double, 3> vecA, std::array<double, 3> vecB);

std::array<int, 3> operator+ (std::array<int, 3> vecA, std::array<int, 3> vecB);

std::array<double, 3> operator- (std::array<double, 3> vecA, std::array<double, 3> vecB);

std::array<double, 3> operator- (std::array<double, 3> vecA, std::array<int, 3> vecB);

std::array<double, 3> operator*(double a, const std::array<double, 3>& vec);

std::array<double, 3> operator*(const std::array<double, 3>& vec, double a);

std::array<double, 3> operator*(const std::array<int, 3>& vec, double a);

double operator* (const std::array<double, 3>& vecA, const std::array<double, 3>& vecB);

std::array<double, 3> operator/ (const std::array<double, 3>& vec, double a);

double abs (const std::array<double, 3> vec);

#endif
