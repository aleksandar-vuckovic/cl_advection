#ifndef SRC_STREAMLINES_HPP_
#define SRC_STREAMLINES_HPP_

#include<array>
#include "Field.hpp"
#include "VelocityField.hpp"


class Streamlines : Field<double> {
private:
	array<vector<array<double, 2>>, 10> streamlines;
public:
	Streamlines(Field<double>& field, VelocityField& vel, double dt);
	virtual ~Streamlines();
};

#endif /* SRC_STREAMLINES_HPP_ */
