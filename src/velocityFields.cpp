/**
 * A library for the definitions of the various velocity fields and its gradients.
 */

#include "velocityFields.hpp"
#include "vecMath.hpp"

/**
 * The functional definition of the shear field
 * It evalutes the shear field at the given coordinates with the given parameter.
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The scaling factor of the shear field
 * @return The velocity vector at the given point
 */
Vector shearField(double x, double y, double z, double v0) {
    Vector tempReturn = {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y), 0};
	return v0*tempReturn;
}

/**
 * The functional defintion of the jacobian matrix of the shear field
 * It evaluates the jacobian matrix of the shear field a tthe given coordinates with the given parameter
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The scaling factor of the shear field
 * @return The jacobian matrix at the given point
 */
Matrix gradShearField(double x, double y, double z, double v0) {
    Matrix tempReturn;
    tempReturn[0] = {-M_PI*cos(M_PI*x)*cos(M_PI*y), M_PI*sin(M_PI*x)*sin(M_PI*y), 0};
    tempReturn[1] = {-M_PI*sin(M_PI*x)*sin(M_PI*y), M_PI*cos(M_PI*x)*cos(M_PI*y), 0};
    tempReturn[2] = {0, 0, 0};

    tempReturn[0] = v0*tempReturn[0];
    tempReturn[1] = v0*tempReturn[1];
    return tempReturn;
}

/**
 * The functional definition of the navier field
 * It evalutes the navier field at the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The value of the x-component at the origin
 * @param c1 A parameter of the field
 * @param c2 A parameter of the field
 * @return The velocity vector at the given point
 */
Vector navierField(double x, double y, double z, double v0, double c1, double c2) {
	return { v0 + c1*x + c2*y, -c1*y, 0};
}

/**
 * The functional definition of the jacobian matrix of the navier field
 * It evalutes the navier field at the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The value of the x-component at the origin
 * @param c1, c2 Parameters of the navier field
 * @return The jacobian matrix at the given point
 */
Matrix gradNavierField(double x, double y, double z, double v0, double c1, double c2) {
    Matrix tempReturn;
	tempReturn[0] = {c1, c2, 0};
	tempReturn[1] = {0, -c1, 0};
	tempReturn[2] = {0, 0, 0};
	return tempReturn;
}

/**
 * The functional definition of the quadratic field.
 * It evaluates the quadratic field a the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param c1 The value of the x-component at the origin
 * @param c2, c3, c4 Parameters of the quadratic of the field
 * @return The velocity vector at the given point
 */
Vector quadraticField(double x, double y, double z, double v0, double c1, double c2, double c3) {
    return {v0 + c1*x + c2*y + c3*y*y, -c1*y};
}

/**
 * The functional definition of the jacobian matrix of the quadratic field
 * It evalutes the navier field at the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param c1 The value of the x-component at the origin
 * @param c2, c3, c4 Parameters of the quadratic field
 * @return The jacobian matrix at the given point
 */
Matrix gradQuadraticField(double x, double y, double z, double v0, double c1, double c2, double c3) {
    Matrix tempReturn;
    tempReturn[0] = {c1, c2 + 2*c3*y, 0};
    tempReturn[1] = {0, -c1, 0};
    tempReturn[2] = {0, 0, 0};
    return tempReturn;
}

/**
 * The functional definition of the strawberry field, a 3D linear field.
 * It evaluates the quadratic field a the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param c1 The value of the x-component at the origin
 * @param c2, c3, c4 Parameters of the quadratic of the field
 * @return The velocity vector at the given point
 */
Vector strawberryField(double x, double y, double z, double v0, double w0, double x0, double y0, double z0,
                     double c1, double c2, double c3, double c4, double c5, double c6) {
    return {v0 + c1 * (x - x0) + c2*y + c3*(z - z0), -(c1 + c6)*y, w0 + c4*(x - x0) + c5*(y - y0) + c6*(z - z0) };
}

Matrix gradStrawberryField(double x, double y, double z, double v0, double w0, double x0, double y0, double z0,
                           double c1, double c2, double c3, double c4, double c5, double c6) {
    Matrix tempReturn;
    tempReturn[0] = {c1, c2, c3};
    tempReturn[1] = {0, -(c1 + c6), 0};
    tempReturn[2] = {c4, c5, c6};
    return tempReturn;
}

