#include "vecMath3D.hpp"
#include "fields.hpp"
#include "BoundaryCondition.hpp"

std::array<double, 3> shearField(double x, double y, double z) {
	return {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y), 0};

}

std::array<std::array<double, 3>, 3> gradShearField(double x, double y, double z) {
    std::array<std::array<double, 3>, 3> tempReturn;
    tempReturn[0] = {-M_PI*cos(M_PI*x)*cos(M_PI*y), -M_PI*sin(M_PI*x)*sin(M_PI*y), 0};
    tempReturn[1] = {-M_PI*sin(M_PI*x)*sin(M_PI*y), M_PI*cos(M_PI*x)*cos(M_PI*y), 0};
    tempReturn[2] = {0, 0, 0};
    return tempReturn;
}
