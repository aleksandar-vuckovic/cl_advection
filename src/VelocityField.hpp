#ifndef CLASS_VELOCITYFIELD
#define CLASS_VELOCITYFIELD
#include <array>

using std::array;

class VelocityField {
private:
	double v0, c1, c2, tau;
	std::string name;
	double maxNormValue;
public:
	VelocityField(std::string name, double v0, double c1, double c2, double tau,
			double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double dx);
	array<double, 3> at(double t, double x, double y, double z);
	array<array<double, 3>, 3> gradAt(double t, double x, double y, double z);
	double getC1();
	double getTau();
	std::string getName();
	double getMaxNormValue();
};

#endif
