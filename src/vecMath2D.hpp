#ifndef VEC_MATH_2D
#define VEC_MATH_2D

#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

std::array<double, 2> operator+ (std::array<double, 2> vecA, std::array<double, 2> vecB);

std::array<double, 2> operator- (std::array<double, 2> vecA, std::array<double, 2> vecB);

std::array<double, 2> operator- (std::array<double, 2> vecA, std::array<int, 2> vecB);

std::array<double, 2> operator*(double a, const std::array<double, 2>& vec);

std::array<double, 2> operator*(const std::array<double, 2>& vec, double a);

std::array<double, 2> operator*(const std::array<int, 2>& vec, double a);

double operator* (const std::array<double, 2>& vecA, const std::array<double, 2>& vecB);

double abs (const std::array<double, 2> vec);

#endif
