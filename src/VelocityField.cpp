#include <string>
#include <array>
#include "VelocityField.hpp"
#include "fields.hpp"

using std::array;

VelocityField::VelocityField(std::string name, double v0, double c1, double c2) {
	this->v0 = v0;
	this->c1 = c1;
	this->c2 = c2;

	if (name == "shearField") {
		field = "shearField";
	} else if (name == "navierField") {
		field = "navierField";
	} else {
		throw std::invalid_argument("No available field was chosen.");
	}
}

array<double, 3> VelocityField::at(double x, double y, double z) {
	if (field == "shearField") {
		return shearField(x, y, z, v0);
	} else if (field == "navierField") {
		return navierField(x, y, z, v0, c1, c2);
	}
	return {0, 0, 0};
}

array<array<double, 3>, 3> VelocityField::gradAt(double x, double y, double z) {
	if (field == "shearField") {
		return gradShearField(x, y, z, v0);
	} else if (field == "navierField") {
		return gradNavierField(x, y, z, v0, c1, c2);
	}
	return {0, 0, 0};
}
