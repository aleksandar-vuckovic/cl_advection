#include <string>
#include <array>
#include "VelocityField.hpp"
#include "velocityFields.hpp"


#define _USE_MATH_DEFINES
#include <cmath>


using std::array;

VelocityField::VelocityField(std::string name, double v0, double c1, double c2, double tau,
		double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double dx, double dy, double dz) {
	this->v0 = v0;
	this->c1 = c1;
	this->c2 = c2;
	this->tau = tau;

	this->xmin = xmin;
	this->xmax = xmax;
	this->ymin = ymin;
	this->ymax = ymax;
	this->zmin = zmin;
	this->zmax = zmax;
	this->dx = dx;
	this->dy = dy;
	this->dz = dz;

	double maxNormValue = 0, currentVal = 0, x, y, z;

	for (int i = 0; i < (xmax - xmin)/dx; i++) {
		for (int j = 0; j < (ymax - ymin)/dx; j++) {
			for (int k = 0; k < (zmax - zmin)/dx; k++) {
				x = i*dx - xmin;
				y = j*dy - ymin;
				z = k*dz - zmin;
				if (name == "shearField")
					currentVal = abs(shearField(x, y, z, v0));
				else if (name == "navierField")
					currentVal = abs(navierField(x, y, z, v0, c1, c2));
				else if (name == "timeDependentNavierField") {
					// Since |cos(x)| <= 1 we assume the worst case and take the maximum absolute value for the norm value of the velocity field
					currentVal = abs(navierField(x, y, z, v0, c1, c2));
				}
				if (currentVal > maxNormValue)
					maxNormValue = currentVal;
			}
		}
	}

	this->maxNormValue = maxNormValue;

	if (name == "shearField") {
		this->name = "shearField";
	} else if (name == "navierField") {
		this->name = "navierField";
	} else if (name == "timeDependentNavierField") {
		this->name = "timeDependentNavierField";
	} else {
		throw std::invalid_argument("No available field was chosen.");
	}
}

array<double, 3> VelocityField::at(double t, double x, double y, double z) {
	if (name == "shearField") {
		return shearField(x, y, z, v0);
	} else if (name == "navierField") {
		return navierField(x, y, z, v0, c1, c2);
	} else if (name == "timeDependentNavierField") {
		return cos(M_PI*t/tau)*navierField(x, y, z, v0, c1, c2);
	}
	return {0, 0, 0};
}

array<array<double, 3>, 3> VelocityField::gradAt(double t, double x, double y, double z) {
	if (name == "shearField") {
		return gradShearField(x, y, z, v0);
	} else if (name == "navierField") {
		return gradNavierField(x, y, z, v0, c1, c2);
	} else if (name == "timeDependentNavierField") {
		return cos(M_PI*t/tau)*gradNavierField(x, y, z, v0, c1, c2);
	}
	return {0, 0, 0};
}

void VelocityField::writeToFile(double t) {
	int numX = (xmax - xmin)/dx;
	int numY = (ymax - ymin)/dy;
	int numZ = (zmax - zmin)/dz;

	double *fieldValues = new double[numX*numY*numZ*3];
	double x, y, z;
	int index = 0;

	for (int k = 0; k < numZ; k++) {
		for (int j = 0; j < numY; j++) {
			for (int i = 0; i < numX; i++) {
				x = i*dx - xmin;
				y = j*dy - ymin;
				z = k*dz - zmin;
				array<double, 3> temp = this->at(t, x, y, z);
				fieldValues[index] = temp[0];
				fieldValues[index + 1] = temp[1];
				fieldValues[index + 2] = temp[2];
				index += 3;
			}
		}
	}

	std::string filename = "data/Vel_t=" + std::to_string(t) + ".bin";
	FILE *velFile;
	velFile = fopen(filename.data(), "wb");
	fwrite(fieldValues, sizeof(double), numX*numY*numZ*3, velFile);
	fclose(velFile);

	delete[] fieldValues;
}

double VelocityField::getC1() {
	return c1;
}

double VelocityField::getTau() {
	return tau;
}

std::string VelocityField::getName() {
	return name;
}

double VelocityField::getMaxNormValue() {
	return maxNormValue;
}
