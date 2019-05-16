#ifndef CLASS_LEVELSET
#define CLASS_LEVELSET

#include <sstream>    // Stringstream
#include <fstream>    // Filestream

#include "Field.hpp"
#include "VelocityField.hpp"
#include "vecMath3D.hpp"

class LevelSet : Field<double> {
private:
    double dx;
    double dy;
    double dz;
    VelocityField *field;

public:
    LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field) : Field<double>(numX, numY, numZ) {
		this->dx = dx;
		this->dy = dy;
		this->dz = dz;
		this->field = field;
    }

    array<double, 3> getInitCP(double dt, array<double, 3> expcp, double epsilon);
    array<double, 3> getContactPoint(double dt, int timestep, int timesteps, array<double, 3> initCP);
    array<int, 3> getContactPointCoordinates(array<double, 3> point);
    double getReferenceAngleExplicitEuler(double dt, int timestep, array<double, 3> n_sigma_init, array<double, 3> CP);
    double getReferenceAngleLinearField(double t, double c1, double c2, double theta0);
    double getContactAngle(double dt, double timestep, array<int, 3> cell);
    double getReferenceCurvature(double dt, double timestep, double initCurvature, array<double, 3> CP, array<int, 3> cell);
    //For now, getCurvature() only works for the stationary droplet
    double getCurvature(double dt, int timestep,  array<int, 3> cell) const;
    double sumLevelSet();
    void writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile);
    void initDroplet(array<double, 3> center, double radius, double epsilon);
    void calculateNextTimestep(double dt, int timestep);
};

#endif
