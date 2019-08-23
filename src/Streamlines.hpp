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
public:
	Streamlines(int numX, int numY, int numZ, VelocityField& vel, double dt);
	void writeToFile();
};

#endif /* SRC_STREAMLINES_HPP_ */
