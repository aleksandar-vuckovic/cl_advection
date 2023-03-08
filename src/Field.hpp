/**
 * @class Field
 * The template for the abstract field class.
 *
 * Save the attribute T for each point on the grid
 */
#ifndef CLASS_FIELD
#define CLASS_FIELD
#include <vector>

template <class T>
class Field
{
private:
    // Used to store the data of the field
    std::vector<T> data;

protected:
    /// The number of cells in each direction
    ///@{
    int numX, numY, numZ;
    ///@}

    double dx, dy, dz;

public:
    Field(int numX, int numY, int numZ, double dx, double dy, double dz);
    T &at(int i, int j, int k);
    const T &at(int i, int j, int k) const;
    const std::vector<T> &getData() const;
};

#include "Field.cpp"

#endif
