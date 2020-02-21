#ifndef FIELD_CPP
#define FIELD_CPP
#include <cmath>
#include "Field.hpp"

/**
 * The constructor.
 *
 * Instantiated the one-dimensional vector with the needed length and sets
 * the member attributes numX, numY and numZ.
 *
 * @param numX, numY, numZ The number of cells in each direction
 */
template<class T>
Field<T>::Field(int numX, int numY, int numZ, double dx, double dy, double dz) : data(numX*numY*numZ){
    this->numX = numX;
    this->numY = numY;
    this->numZ = numZ;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
}

/**
 * Evalutes the field at the given indices
 *
 * @param x, y, z The indices of the point
 * @return The value of T at the point
 */
template<class T>
T& Field<T>::at(int i, int j, int k) {
    return data[i + j*numX + k*numX*numY];
}

/**
 * Evalutes the field at the given indices.
 * This function is called when the field is evaluated within another function defined as "const".
 * Further, it is used when a const variable is set to the value of a point of the field.
 *
 * @param x, y, z The indices of the point
 * @return The value of T a the point
 */
template<class T>
const T& Field<T>::at(int i, int j, int k) const {
    return data[i + j*numX + k*numX*numY];
}

template<class T>
const std::vector<T>& Field<T>::getData() const {
	return data;
}

#endif /* FIELD_CPP */
