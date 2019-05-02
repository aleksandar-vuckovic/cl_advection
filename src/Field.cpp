#ifndef FIELD_CPP
#define FIELD_CPP
#include "Field.hpp"
template<class T>
Field<T>::Field(int numX, int numY, int numZ) {
    data = std::vector<T>(numX*numY*numZ);
    this->numX = numX;
    this->numY = numY;
    this->numZ = numZ;
}

template<class T>
T& Field<T>::at(int x, int y, int z) {
    return data[x + y*numX + z*numX*numY];
}

template<class T>
const T& Field<T>::at(int x, int y, int z) const {
    return data[x + y*numX + z*numX*numY];
}

#endif /* FIELD_CPP */
