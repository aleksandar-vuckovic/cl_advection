/**
 * A library for the definitions of the various velocity fields and its gradients.
 */

#include "vecMath3D.hpp"
#include "velocityFields.hpp"

std::array<double, 3> shearField(double x, double y, double z, double v0) {
	std::array<double, 3> tempReturn = {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y), 0};
	return v0*tempReturn;
}

std::array<std::array<double, 3>, 3> gradShearField(double x, double y, double z, double v0) {
    std::array<std::array<double, 3>, 3> tempReturn;
    tempReturn[0] = {-M_PI*cos(M_PI*x)*cos(M_PI*y), M_PI*sin(M_PI*x)*sin(M_PI*y), 0};
    tempReturn[1] = {-M_PI*sin(M_PI*x)*sin(M_PI*y), M_PI*cos(M_PI*x)*cos(M_PI*y), 0};
    tempReturn[2] = {0, 0, 0};

    tempReturn[0] = v0*tempReturn[0];
    tempReturn[1] = v0*tempReturn[1];
    return tempReturn;
}

std::array<double, 3> navierField(double x, double y, double z, double v0, double c1, double c2) {
	return { v0 + c1*x + c2*y, -c1*y, 0};
}

std::array<std::array<double, 3>, 3> gradNavierField(double x, double y, double z, double v0, double c1, double c2) {
	std::array<std::array<double, 3>, 3> tempReturn;
	tempReturn[0] = {c1, c2, 0};
	tempReturn[1] = {0, -c1, 0};
	tempReturn[2] = {0, 0, 0};
	return tempReturn;
}
