#ifndef CLASS_FIELD
#define CLASS_FIELD
#include <vector>

template <class T>
class Field {
    private:
    std::vector<T> data;

protected:
    int numX, numY, numZ;

public:
    Field(int numX, int numY, int numZ);
    T& at(int x, int y, int z);
    const T& at(int x, int y, int z) const;

};
#include "Field.cpp"

#endif
