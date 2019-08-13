#ifndef CLASS_LEVELSET
#define CLASS_LEVELSET

#include <sstream>    // Stringstream
#include <fstream>    // Filestream

#include "Field.hpp"
#include "VelocityField.hpp"
#include "vecMath3D.hpp"

class LevelSet : Field<double> {
private:
    //The width of a cell in each direction
    double dx, dy, dz;

    // A pointer to the VelocityField acting on the LevelSet field.
    VelocityField *field;

    // Decides which contact point to track. Only applicable in 2D.
    std::string trackedCP;

    /**
    *   A pointer to reference data calculated with the explicit euler algorithm.
    *   This is needed, since many other reference solvers are coupled, and all of them require this data,
    *   leading to a high number in nested loops.
    **/
    std::vector< array<double, 3> > *positionReference;

public:
    LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field, std::string trackedCP, std::vector< array<double, 3> > *positionReference);

    array<double, 3> getInitCP(array<double, 3> expcp, double epsilon);
    array<double, 3> getContactPointExplicitEuler(double dt, int timestep, array<double, 3> initCP);
    array<double, 3> getContactPointLinearField(double t, double c1, double x0, double v0);
    array<int, 3> getContactPointIndices(array<double, 3> point);
    array<double, 3> getContactPoint(array<int, 3> indices);
    double getReferenceAngleExplicitEuler(double dt, int timestep, array<double, 3> n_sigma_init, array<double, 3> CP_init);
    double getReferenceAngleLinearField(double t, double c1, double c2, double theta0);
    double getContactAngle(array<int, 3> cell);
    double getReferenceCurvatureExplicitEuler(double dt, int timestep, double initCurvature, double initAngle, array<double, 3> CP);
    double getReferenceCurvatureLinearField(double t, double init_curvature);
    double getReferenceCurvatureQuadraticField(double t, double init_curvature);
    double getCurvatureDivergence(array<int, 3> cell) const;
    double getCurvatureHeight(array<int, 3> cell) const;
    void writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile);
    void writeNormalVectorToFile(double t);
    double sumLevelSet();
    void initDroplet(array<double, 3> center, double radius);
    void calculateNextTimestep(double dt, int timestep);
};

#endif
