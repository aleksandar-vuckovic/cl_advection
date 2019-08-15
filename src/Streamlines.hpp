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
	std::array<std::vector<std::array<double, 3>>, 10> streamlines;
public:
	Streamlines(int numX, int numY, int numZ, VelocityField& vel, double dt);
	void writeToFile();
	virtual ~Streamlines();
};

#endif /* SRC_STREAMLINES_HPP_ */
