#ifndef CLASS_LEVELSET
#define CLASS_LEVELSET

#include <sstream>    // Stringstream
#include <fstream>    // Filestream

#include "Field.hpp"
#include "vecMath.hpp"
#include "VelocityField.hpp"
#include "enums.hpp"

#include <omp.h>

class LevelSet : Field<double> {
private:
    // A pointer to the VelocityField acting on the LevelSet field.
    VelocityField* field;

    // Decides which contact point to track. Only applicable in 2D.
    std::string trackedCP, outputDirectory;

    /**
    *   Reference data..
    *   This is needed, since many other reference solvers are coupled, and all of them require this data,
    *   leading to a high number in nested loops. 
    **/
    std::vector<Vector> positionReference;
    std::vector<Vector> normalReference;
    std::vector<double> angleReference;
    std::vector<double> curvatureReference;

public:
    LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field,
            std::string trackedCP, double dt, int timesteps, Vector expCP, Vector expNormalVec, double expAngle,
            double initCurvature, std::string outputDirectory);

    Vector getInitCP(Vector expcp, double epsilon);
    void contactPointExplicitEuler(double dt, int timesteps, Vector initCP);
    void contactPointLinearField(double dt, int timesteps, double c1, double x0, double v0);
    Vector getContactPoint(int timestep, bool indexOnly = false) const;
    array<int, 3> getContactPointIndices(int timestep) const;
    void referenceNormalExplicitEuler(double dt, int timestep, Vector n_sigma_init);
    void referenceAngleLinearField(double dt, int timesteps, double theta0);
    Vector getNormalVector(array<int, 3> cell,  bool useInterpolation = true) const;
    Vector getNormalVector(int i, int j, int k) const;
    Vector getTangentialVector(Vector normal) const;
    double getContactAngleInterpolated(int timestep);
    void referenceCurvatureExplicitEuler2D(double dt, int timestep, double initCurvature);
    double referenceCurvatureDeriv3D(double initCurvature, Matrix expNormalVectorGradient);
    void referenceCurvatureLinearField(double dt, int timesteps, double initCurvature);
    void referenceCurvatureQuadraticField(double dt, int timesteps, double initCurvature);
    double getCurvature(array<int, 3> cell) const;
    double getCurvature(int i, int j, int k) const;
    double getCurvatureInterpolated(int timestep) const;
    void writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *mainXmfFile, std::ofstream *tauXmfFile);
    void writeTangentialVectorToFile(double dt, int timestep);
    double sumLevelSet();
    void initSphere(Vector center, double radius);
    void initPlane(Vector refPoint, double polarAngle, double azimuthalAngle);
    void initParaboloid(Vector refPoint, double stretchX, double stretchY, double heightMinimum);
    void initEllipsoid(Vector refPoint, double stretchX, double stretchY, double stretchZ);
    static Vector expectedNormalVector(Vector contactPoint, InitShape shape, Vector refPoint, std::vector<double> params);
    static Matrix expectedNormalVectorGradient(Vector contactPoint, InitShape shape, Vector refPoint, std::vector<double> params);
    Vector normalVector2D(double initAngle);
    static Vector normalVector2D(double initAngle, std::string trackedCP);
    void calculateNextTimestep(double dt, int timestep);

    std::vector<Vector> getPositionReference() const;
    std::vector<double> getAngleReference() const;
    const std::vector<double>& getCurvatureReference() const;
};

Vector normalVector2D(double initAngle, std::string trackedCP);

#endif
