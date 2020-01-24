#ifndef CLASS_VELOCITYFIELD
#define CLASS_VELOCITYFIELD

#include <string>
#include <array>
#include "velocityFields.hpp"

using std::array;

class VelocityField {
private:
    //@{
    /** The space the velocity field is defined on */
	double xmin, xmax, ymin, ymax, zmin, zmax;
	//@}

	//@{
	/** The width of a cell in each direction */
	double dx, dy, dz;
	//@}

	//@{
	/** The parameters of the velocity field */
    double v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6, tau, azimuthalAngle;
    Matrix rotMatrix;
	//@}

    /// The kind of velocity field, either navier, navier with cosine modulation, shear or strawberry.
    std::string name, outputDirectory;
	///The maximum absolute value of the field on the space it is defined on.
	double maxAbsoluteValue;
public:
    VelocityField(std::string name, double v0, double w0, double x0, double y0, double z0,
                                 double c1, double c2, double c3, double c4, double c5, double c6, double tau,
                                 double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                                 double dx, double dy, double dz, double azimuthalAngle, std::string outputDirectory);
    Vector at(double t, double x, double y, double z);
    Matrix gradAt(double t, double x, double y, double z);
	void writeToFile(double t);
	double getXMax();
	double getYMax();
	double getDx();
	double getDy();
	double getV0();
	double getC1();
	double getC2();
	double getC3();
	double getTau();
	std::string getName();
	double getMaxNormValue();
    double getAzimuthalAngle();
};

#endif
