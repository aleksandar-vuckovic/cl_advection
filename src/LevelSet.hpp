#ifndef CLASS_LEVELSET
#define CLASS_LEVELSET

#include <sstream>    // Stringstream
#include <fstream>    // Filestream

#include "Field.hpp"
#include "vecMath.hpp"
#include "VelocityField.hpp"

#include <omp.h>

class LevelSet : Field<double> {
private:
    //The width of a cell in each direction
    double dx, dy, dz, dt;

    // A pointer to the VelocityField acting on the LevelSet field.
    VelocityField *field;

    // Decides which contact point to track. Only applicable in 2D.
    std::string trackedCP;

    /**
    *   Reference data..
    *   This is needed, since many other reference solvers are coupled, and all of them require this data,
    *   leading to a high number in nested loops.
    **/
    std::vector< array<double, 3>> positionReference;
    std::vector<double> angleReference;
    std::vector<double> curvatureReference;

public:
    LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field, std::string trackedCP, double dt, int timesteps,
            array<double, 3> expcp, double expAngle, double initCurvature);

    array<double, 3> getInitCP(array<double, 3> expcp, double epsilon);
    void contactPointExplicitEuler(double dt, int timestep, array<double, 3> initCP);
    array<double, 3> getContactPointLinearField(double t, double c1, double x0, double v0);
    array<int, 3> getContactPointIndices(array<double, 3> point);
    array<double, 3> getContactPoint(array<int, 3> indices);
    void referenceAngleExplicitEuler(double dt, int timestep, array<double, 3> n_sigma_init, array<double, 3> CP_init);
    void referenceAngleLinearField(double dt, int last_timestep, double theta0);
    double getContactAngle(array<int, 3> cell);
    void referenceCurvatureExplicitEuler(double dt, int timestep, double initCurvature, double initAngle, array<double, 3> CP);
    void referenceCurvatureLinearField(double dt, int timesteps, double init_curvature);
    void referenceCurvatureQuadraticField(double dt, int timesteps, double init_curvature);
    double getCurvatureDivergence(array<int, 3> cell) const;
    double getCurvatureHeight(array<int, 3> cell) const;
    void writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile);
    void writeTangentialVectorToFile(double t);
    double sumLevelSet();
    void initDroplet(array<double, 3> center, double radius);
    array<double, 3> normalVector2D(double initAngle);
    void calculateNextTimestep(double dt, int timestep);

    std::vector<array<double, 3>> getPositionReference();
    std::vector<double> getAngleReference();
    std::vector<double> getCurvatureReference();
};

#endif
