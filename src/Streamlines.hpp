#ifndef SRC_STREAMLINES_HPP_
#define SRC_STREAMLINES_HPP_

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include "LevelSet.hpp"
#include "VelocityField.hpp"
#include "vecMath.hpp"


class Streamlines : Field<int> {
private:
    std::string outputDirectory;
public:
    Streamlines(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField* field, double dt, std::string outputDirectory);
	void writeToFile();
};

#endif /* SRC_STREAMLINES_HPP_ */
