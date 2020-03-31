#ifndef FIELDS_HPP_
#define FIELDS_HPP_
#include "vecMath.hpp"

Vector shearField(double x, double y, double z, double v0);
Matrix gradShearField(double x, double y, double z, double v0);
array<double, 6> partialsShearField(double x, double y, double z, double v0);
Vector navierField(double x, double y, double z, double v0, double c1, double c2);
Matrix gradNavierField(double x, double y, double z, double v0, double c1, double c2);
Vector quadraticField(double x, double y, double z, double v0, double c1, double c2, double c3);
Matrix gradQuadraticField(double x, double y, double z, double v0, double c1, double c2, double c3);
array<double, 6> partialsQuadraticField(double x, double y, double z, double v0, double c1, double c2, double c3);
Vector strawberryField(double x, double y, double z, double v0, double w0, double x0, double y0, double z0,
                     double c1, double c2, double c3, double c4, double c5, double c6);
Matrix gradStrawberryField(double x, double y, double z, double v0, double w0, double x0, double y0, double z0,
                     double c1, double c2, double c3, double c4, double c5, double c6);

#endif /* FIELDS_HPP_ */
