#ifndef CLASS_VELOCITYFIELD
#define CLASS_VELOCITYFIELD
#include <array>

using std::array;

class VelocityField {
private:
	double v0, c1, c2;
	std::string field;
public:
	VelocityField(std::string name, double v0, double c1, double c2);
	array<double, 3> at(double x, double y, double z);
	array<array<double, 3>, 3> gradAt(double x, double y, double z);
};

#endif
