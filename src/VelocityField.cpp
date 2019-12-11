/**
 * @class VelocityField
 * The class of the velocity field acting on the Level set field.
 *
 * An object of class VelocityField represents a velocity field that knows its type, meaning whether it is
 * a navier or a shear field, its parameters and the space it is defined on.
 */

#include "VelocityField.hpp"

#define _USE_MATH_DEFINES
#include <cmath>


using std::array;

/**
 * The constructor.
 *
 * @param name The kind of the field. This decides the underlying function the velocity field will apply at each point
 * @param v0 For the shear field, this is a scaling factor. For the navier field, this is the velocity of the x-component in the origin.
 * @param \f$c_i\f$ Parameters of the velocity field, currently only used for the navier field
 * @param tau tau/2 is the period of osciallation of the time dependent navier field.
 * @param xmin, xmax, ymin, ymax, zmin, zmax The space the velocity field is defined on
 * @param dx, dy, dz The width of a cell in each direction
 */
VelocityField::VelocityField(std::string name, double v0, double w0, double x0, double y0, double z0,
                             double c1, double c2, double c3, double c4, double c5, double c6, double tau,
                             double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                             double dx, double dy, double dz, double azimuthalAngle) {
	this->v0 = v0;
    this->w0 = w0;
    this->x0 = x0;
    this->y0 = y0;
    this->z0 = z0;
	this->c1 = c1;
	this->c2 = c2;
	this->c3 = c3;
    this->c4 = c4;
    this->c5 = c5;
    this->c6 = c6;
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
    this->azimuthalAngle = azimuthalAngle/180 * M_PI;

    // Matrix for rotation around y-axis
    array<double, 3> row1 = { cos(azimuthalAngle), 0, -sin(azimuthalAngle)};
    array<double, 3> row2 = {0, 1, 0};
    array<double, 3> row3 = {sin(azimuthalAngle), 0, cos(azimuthalAngle)};

    rotMatrix = {row1, row2, row3};

	double maxNormValue = 0, currentVal = 0, x, y, z;

	for (int i = 0; i < (xmax - xmin)/dx; i++) {
		for (int j = 0; j < (ymax - ymin)/dy; j++) {
			for (int k = 0; k < (zmax - zmin)/dz; k++) {
				x = i*dx - xmin;
				y = j*dy - ymin;
				z = k*dz - zmin;
				if (name == "shearField")
                    currentVal = abs(rotMatrix * shearField(x, y, z, v0));
				else if (name == "navierField")
                    currentVal = abs(rotMatrix * navierField(x, y, z, v0, c1, c2));
				else if (name == "timeDependentNavierField") {
					// Since |cos(x)| <= 1 we assume the worst case and take the maximum absolute value for the norm value of the velocity field
                    currentVal = abs(rotMatrix * navierField(x, y, z, v0, c1, c2));
				}
				else if (name == "quadraticField") {
                    currentVal = abs(rotMatrix * quadraticField(x, y, z, v0, c1, c2, c3));
				}
                else if (name == "strawberryField") {
                    currentVal = abs(rotMatrix * strawberryField(x, y, z, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6));
                }
				if (currentVal > maxNormValue)
					maxNormValue = currentVal;
			}
		}
	}

	this->maxAbsoluteValue = maxNormValue;

    if (name == "shearField" || name == "navierField" || name == "timeDependentNavierField" || name == "quadraticField" || name == "strawberryField") {
		this->name = name;
	} else {
		throw std::invalid_argument("No available field was chosen.");
	}
}

/**
 * Evaluates the velocity field at the given point.
 *
 * @param t The time
 * @param x, y, z The coordinates of the point
 * @return The velocity field at the given coordinates
 */
array<double, 3> VelocityField::at(double t, double x, double y, double z) {

	if (name == "shearField") {
		return rotMatrix *  shearField(x, y, z, v0);
	} else if (name == "navierField") {
		return rotMatrix * navierField(x, y, z, v0, c1, c2);
	} else if (name == "timeDependentNavierField") {
		return cos(M_PI*t/tau) * rotMatrix * navierField(x, y, z, v0, c1, c2);
	} else if (name == "quadraticField") {
	    return rotMatrix * quadraticField(x, y, z, v0, c1, c2, c3);
    } else if (name == "strawberryField") {
        return strawberryField(x, y, z, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6);
    }
	return {0, 0, 0};
}

/**
 * Evaluates the jacobian matrix at a given point
 *
 * @param t The time
 * @param x, y, z The coordinates of the point
 * @return The jacobian matrix at the given coordinates
 */
Matrix VelocityField::gradAt(double t, double x, double y, double z) {

	if (name == "shearField") {
        return rotMatrix*gradShearField(x, y, z, v0);
	} else if (name == "navierField") {
        return rotMatrix*gradNavierField(x, y, z, v0, c1, c2);
	} else if (name == "timeDependentNavierField") {
        return cos(M_PI*t/tau)*rotMatrix*gradNavierField(x, y, z, v0, c1, c2);
	} else if (name == "quadraticField") {
        return rotMatrix*gradQuadraticField(x, y, z, v0, c1, c2, c3);
    } else if (name == "strawberryField") {
        return rotMatrix*gradStrawberryField(x, y, z, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6);
    }
	return {0, 0, 0};
}

/**
 * Writes the velocity field at a given time to disk
 *
 * Writes the velocity field for visualization in Paraview. Since no XMF file is written,
 * simply calling this function is not enough for visualization. Thus, this function is called within
 * LevelSet::writeToFile, which does write a XMF file.
 *
 * @param t The time
 */
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

double VelocityField::getXMax() {
	return xmax;
}
double VelocityField::getYMax() {
	return ymax;
}

double VelocityField::getDx() {
	return dx;
}

double VelocityField::getDy() {
	return dy;
}
double VelocityField::getV0() {
    return v0;
}

double VelocityField::getC1() {
	return c1;
}

double VelocityField::getC2() {
    return c2;
}

double VelocityField::getC3() {
    return c3;
}

double VelocityField::getTau() {
	return tau;
}

std::string VelocityField::getName() {
	return name;
}

double VelocityField::getMaxNormValue() {
	return maxAbsoluteValue;
}

double VelocityField::getAzimuthalAngle() {
    return azimuthalAngle;
}
