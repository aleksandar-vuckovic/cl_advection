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
 * @param azimuthalAngle, alpha Parameters defining the rotation of the velocity field. A rotation introduces a new coordinate system (\f$(\tilde{x}, y)\f$) with
 * \f$\tilde{x} = (x, 0, z) \cdot \hat{n}_\Gamma. \hat{n}_\Gamma and \hat{n}_y are the unit vectors of the new coordinate system.
 */
VelocityField::VelocityField(std::string name, double v0, double w0, double x0, double y0, double z0,
                             double c1, double c2, double c3, double c4, double c5, double c6, double tau,
                             double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                             double dx, double dy, double dz, double azimuthalAngle, double alpha, std::string outputDirectory)
                             : alpha(alpha), n_gamma({cos(azimuthalAngle), 0, sin(azimuthalAngle)}), n_y({0, 1, 0})
{
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
    this->outputDirectory = outputDirectory;

    if (name == "shearField" || name == "navierField" || name == "timeDependentNavierField" || name == "quadraticField" || name == "strawberryField") {
		this->name = name;
	} else {
		throw std::invalid_argument("No available field was chosen.");
	}

	double maxNormValue = 0, currentVal = 0, x, y, z;

	for (int i = 0; i < (xmax - xmin)/dx; i++) {
		for (int j = 0; j < (ymax - ymin)/dy; j++) {
			for (int k = 0; k < (zmax - zmin)/dz; k++) {
				x = i*dx - xmin;
				y = j*dy - ymin;
				z = k*dz - zmin;
				
				currentVal = abs(this->at(0, x, y, z));
				if (currentVal > maxNormValue)
					maxNormValue = currentVal;
			}
		}
	}

	this->maxAbsoluteValue = maxNormValue;
}

/**
 * Evaluates the velocity field at the given point.
 *
 * @param t The time
 * @param x, y, z The coordinates of the point
 * @return The velocity field at the given coordinates
 */
Vector VelocityField::at(double t, double x, double y, double z) {

	double x_tilde = x*n_gamma[0] + z*n_gamma[2] - alpha;
	Vector vec;
	if (name == "shearField") {
		vec = shearField(x_tilde, y, z, v0);
        vec = vec[0]*n_gamma + vec[1]*n_y;
	} else if (name == "navierField") {
		vec = navierField(x_tilde, y, z, v0, c1, c2);
        vec = vec[0]*n_gamma + vec[1]*n_y;
	} else if (name == "timeDependentNavierField") {
		vec = cos(M_PI*t/tau) * navierField(x_tilde, y, z, v0, c1, c2);
        vec = vec[0]*n_gamma + vec[1]*n_y;
	} else if (name == "quadraticField") {
	    vec = quadraticField(x_tilde, y, z, v0, c1, c2, c3);
        vec = vec[0]*n_gamma + vec[1]*n_y;
    } else if (name == "strawberryField") {
        vec = strawberryField(x_tilde, y, z, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6);
    }
    return vec;
}

/**
 * Evaluates the jacobian matrix at a given point
 *
 * @param t The time
 * @param x, y, z The coordinates of the point
 * @return The jacobian matrix at the given coordinates
 */
Matrix VelocityField::gradAt(double t, double x, double y, double z) {
	
	double x_tilde = x*n_gamma[0] + z*n_gamma[2] - alpha;
	Matrix mat;
	if (name == "shearField") {
        mat = gradShearField(x_tilde, y, z, v0);
	} else if (name == "navierField") {
        mat = gradNavierField(x_tilde, y, z, v0, c1, c2);
	} else if (name == "timeDependentNavierField") {
        mat = cos(M_PI*t/tau)*gradNavierField(x_tilde, y, z, v0, c1, c2);
	} else if (name == "quadraticField") {
        mat = gradQuadraticField(x_tilde, y, z, v0, c1, c2, c3);
    } else if (name == "strawberryField") {
        return gradStrawberryField(x_tilde, y, z, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6);
    }
	/* const double c = abs(n_gamma);
	const Vector Q_0 {n_gamma[0]*n_gamma[0]/(c*c), n_gamma[0]/c, n_gamma[0]*n_gamma[2]/(c*c)};
	const Vector Q_1 {n_gamma[0]/c, 1, n_gamma[2]/c};
	const Vector Q_2 {n_gamma[0]*n_gamma[2]/(c*c), n_gamma[2]/c, n_gamma[2]*n_gamma[2]/(c*c)};
	const Matrix Q {Q_0, Q_1, Q_2};
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			mat[i][j] *= Q[i][j]; */

	Matrix temp;
	const double c = abs(n_gamma);
	temp[0][0] = mat[0][0] * n_gamma[0] * n_gamma[0] / (c*c);
	temp[0][1] = mat[0][1] * n_gamma[0] / c;
	temp[0][2] = mat[0][0] * n_gamma[0] * n_gamma[2] / (c*c);
	temp[1][0] = mat[1][0] * n_gamma[0] / c;
	temp[1][1] = mat[1][1];
	temp[1][2] = mat[1][0] * n_gamma[2] / c;
	temp[2][0] = mat[0][0] * n_gamma[0] * n_gamma[2] / (c*c);
	temp[2][1] = mat[0][1] * n_gamma[2] / c;
	temp[2][2] = mat[0][0] * n_gamma[2] * n_gamma[2] / (c*c);

	return temp;
}

/* Matrix  VelocityField::hessianAt(double t, double x, double y, double z) {
	Vector p = {x, y, z};
	p = rotMatrix.transposed*p;
	x = p[0]; 
	y = p[1]; 
	z = p[2];

	if (name == "shearField") {
		Vector row1 = {v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1]) + v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1]),
								v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1]) + v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1]),
								0};
		Vector row2 = {-v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1]) - v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1]),
									-v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1]) -v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1]),
								0};
		Vector row3 = {0, 0, 0};
		Matrix M = {row1, row2, row3};
		return rotMatrix*M;
	} else {
		return {0, 0, 0};
	}
} */

/**
 * Writes the velocity field at a given time to disk
 *
 * Writes the velocity field for visualization in Paraview. Since no XMF file is written,
 * simply calling this function is not enough for visualization. Thus, this function is called within
 * LevelSet::writeToFile, which does write a XMF file.
 *
 * @param t The time
 */

Vector VelocityField::secondPartial(double t, double x, double y, double z, Vector tau) {

	double x_tilde = x*n_gamma[0] + z*n_gamma[2] - alpha;

	// Second order partial derivatives of the velocity fields with the structure
	// {dxdx_vx, dxdy_vx, dydy_vx, dxdx_vy, dxdy_vy, dydy_vy}
	array<double, 6> partials; 

	if (name == "shearField") {
        partials = partialsShearField(x_tilde, y, z, v0);
	} else if (name == "quadraticField") {
        partials = partialsQuadraticField(x_tilde, y, z, v0, c1, c2, c3);
	} else {
		partials = {0, 0, 0, 0, 0, 0};
	}

	Matrix temp;
	double n_x = n_gamma[0];
	double n_z = n_gamma[2];
	double c = abs(n_gamma);

	temp[0][0] = pow(n_x/c, 3) * partials[0] * tau[0] 
                   + pow(n_x/c, 2) * partials[1] * tau[1]
                   + pow(n_x, 2) * n_z / pow(c, 3) * partials[0] * tau[2];

	temp[0][1] = pow(n_x/c, 2) * partials[1] * tau[0] 
	           + n_x/c * partials[2] * tau[1] 
	           + n_x*n_z/pow(c, 2) * partials[1] * tau[2];

	temp[0][2] = pow(n_x, 2)*n_z/pow(c, 3) * partials[0] * tau[0] 
	           + n_x*n_z /pow(c, 3) * partials[1] * tau[1] 
	           + n_x*pow(n_z, 2) / pow(c, 3)* partials[0]*tau[2];

	temp[1][0] = pow(n_x/c, 2) * partials[3] * tau[0] 
	           + n_x/c * partials[4] * tau[1] 
	           + n_x*n_z/pow(c, 2) * partials[3] * tau[2];

	temp[1][1] = n_x/c * partials[4] * tau[0] 
	           + partials[5] * tau[1] 
	           + n_z/c * partials[4] * tau[2];

	temp[1][2] = n_x*n_z/pow(c, 2) *partials[3] * tau[0] 
	           + n_z/c * partials[4] * tau[1] 
	           + pow(n_z/c, 2) * partials[3] * tau[2];

	temp[2][0] = pow(n_x, 2)*n_z/pow(c, 3) * partials[0] * tau[0] 
	           + n_x*n_z/pow(c, 2)*partials[1] * tau[1] 
	           + n_x*pow(n_z, 2)/ pow(c, 3) * partials[0] * tau[2];

	temp[2][1] = n_x*n_z/pow(c, 2) * partials[1] * tau[0] 
	           + n_z/c * partials[2] * tau[1] 
	           + pow(n_z/c, 2) * partials[1] * tau[2];

	temp[2][2] = n_x*pow(n_z, 2)/pow(c, 3) * partials[0] * tau[0] 
	           + pow(n_z/c, 2) * partials[1] * tau[1] 
                   + pow(n_z/c, 3) * partials[0] * tau[2];

	return temp*tau;
	
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
				Vector temp = this->at(t, x, y, z);
				fieldValues[index] = temp[0];
				fieldValues[index + 1] = temp[1];
				fieldValues[index + 2] = temp[2];
				index += 3;
			}
		}
	}

    std::string filename = outputDirectory + "data/Vel_t=" + std::to_string(t) + ".bin";
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
