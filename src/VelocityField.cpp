#include <string>
#include <array>
#include "VelocityField.hpp"
#include "velocityFields.hpp"

using std::array;

VelocityField::VelocityField(std::string name, double v0, double c1, double c2,
		double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double dx) {
	this->v0 = v0;
	this->c1 = c1;
	this->c2 = c2;

	double maxNormValue = 0, currentVal = 0, x, y, z;

	for (int i = 0; i < (xmax - xmin)/dx; i++) {
		for (int j = 0; j < (ymax - ymin)/dx; j++) {
			for (int k = 0; k < (zmax - zmin)/dx; k++) {
				x = i*dx - xmin;
				y = j*dx - ymin;
				z = k*dx - zmin;
				if (name == "shearField")
					currentVal = abs(shearField(x, y, z, v0));
				else if (name == "navierField")
					currentVal = abs(navierField(x, y, z, v0, c1, c2));
				if (currentVal > maxNormValue)
					maxNormValue = currentVal;
			}
		}
	}

	this->maxNormValue = maxNormValue;

	if (name == "shearField") {
		this->field = "shearField";
	} else if (name == "navierField") {
		this->field = "navierField";
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

double VelocityField::getMaxNormValue() {
	return maxNormValue;
}
