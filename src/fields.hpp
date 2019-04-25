#ifndef FIELDS_HPP_
#define FIELDS_HPP_
#include "vecMath3D.hpp"

std::array<double, 3> shearField(double x, double y, double z);
std::array<std::array<double, 3>, 3> gradShearField(double x, double y, double z);


#endif /* FIELDS_HPP_ */
