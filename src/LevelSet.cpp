/**
 * @class LevelSet
 * The class of the Level set field.
 *
 * This class contains most of the methods used in the program, such as calculating the reference contact angle,
 * calculating the contact angle from the simulation data, tracking the contact point and evolving the field with time.
 */

#include "LevelSet.hpp"

/**
 * The constructor.
 * @param numX The number of cells in x direction
 * @param dx Width of one cell in x direction
 * @param field Pointer to the velocity field object acting on the levelset field
 * @param trackedCP Which contact point to track. Only applicable in 2D.
 */
LevelSet::LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field,
        std::string trackedCP, double dt, int timesteps, Vector expCP, double expAngle,
        double initCurvature, InitShape shape, std::vector<double> shapeParams, Vector initCenter, std::string outputDirectory)
        : Field<double>(numX, numY, numZ, dx, dy, dz),
          positionReference(timesteps), normalReference(timesteps),
          tangentAReference(timesteps),tangentBReference(timesteps), tangentCReference(timesteps),
          angleReference(timesteps), curvatureReference(timesteps),
          sectionalCurvatureAReference(timesteps), sectionalCurvatureBReference(timesteps), sectionalCurvatureCReference(timesteps),
          alpha1(timesteps), alpha2(timesteps),
          shape(shape), shapeParams(shapeParams), initCenter(initCenter) {
		this->field = field;
		this->trackedCP = trackedCP;
        this->outputDirectory = outputDirectory;
        expAngle = expAngle/180*M_PI;
        Vector expNormalVec = expectedNormalVector(expCP);

		// Calculate reference data
        if (numZ == 1 && ( field->getName() == "navierField" || field->getName() == "timeDependentNavierField") )
            contactPointLinearField(dt, timesteps, field->getC1(), expCP[0], field->getV0());
        else
            referenceContactPointExplicitEuler(dt, timesteps, expCP);

        referenceNormalExplicitEuler(dt, timesteps, expNormalVec);
        referenceTangentExplicitEuler(dt, timesteps, expNormalVec);
        referenceCurvatureExplicitEuler(dt, timesteps, initCurvature, "sectionalCurvature");
}

/**
 * Calculate the initial contact point from a given expected contact point.
 * Using the expected contact point, iterate over all points of the field, looking for the closest point to expcp which
 * has an absolute Level set value of less than epsilon and lies on the interface.
 *
 * @param expcp The expected contact point
 * @param epsilon Maximum absolute Level set value to consider a point "on the interface"
 * @return The contact point closest to the expected contact point
 *
 */
Vector LevelSet::getInitCP(Vector expcp, double epsilon) {
    Vector candidate = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
        for (int y = 0; y < this->numY; y++)
            for (int z = 0; z < this->numZ; z++) {
                Vector other = {x*dx, y*dy, z*dz};
                if (abs(candidate - expcp) > abs(other - expcp) && std::abs(this->at(x, y, z)) < epsilon && y == 0) {
                    candidate = other;
		}
	}
    return candidate;
}

/**
 * Calculate the new position of the contact point at a given time.
 * Using the explicit euler method, calculate the position of the contact point initCP at all times and writes them to
 * the reference array. If the contact p
 *
 * @param dt The length of a single timestep
 * @param timestep The index of the timestep
 * @param initCP The position of the contact point at timestep 0
 */
void LevelSet::referenceContactPointExplicitEuler(double dt, int last_timestep, Vector initCP) {
    Vector &temp = initCP;
    int factor = 100;
    dt = dt / factor;
    for (int i = 0; i < factor * last_timestep; i++) {
        if (i % factor == 0)
            positionReference[i/factor] = temp;
        temp = temp + dt*field->at(i*dt, temp[0], temp[1], temp[2]);
    }
}

std::vector<Vector> LevelSet::referenceContactPointBackwards(double dt, int last_timestep, Vector point) {
    Vector &temp = point;
    std::vector<Vector> tempReturn;
    int factor = 100;
    dt = dt / factor;
    tempReturn.push_back(temp);
    for (int i = factor * last_timestep; i > 0; i--) {
        temp = temp - dt*field->at((i-1)*dt, temp[0], temp[1], temp[2]);
        tempReturn.push_back(temp);
    }
    return tempReturn;
}

/**
 * Calculate the contact point position for the navier field using the analytic solution.
 *
 * @param t The time at which to calculate the contact point position
 * @param c1 A parameter of the navier field
 * @param x0 The initial position of the contact point
 * @param v0 A parameter of the navier field
 * @return The position of the contact point in arbitrary units
 */
void LevelSet::contactPointLinearField(double dt, int last_timestep, double c1, double x0, double v0) {
    double t = 0;
    for (int i = 0; i < last_timestep; i++) {
        t = i*dt;
        if (field->getName() == "navierField") {
            positionReference[i] = {x0 * exp(c1 * t) + v0/c1 * (exp(c1 * t) - 1), 0, 0};
        } else if (field->getName() == "timeDependentNavierField") {
            double tau = field->getTau();
            positionReference[i] =  {x0 * exp(c1/M_PI * tau*sin(M_PI*t/tau)) + v0/c1 * (exp(c1/M_PI * tau*sin(M_PI*t/tau)) -1), 0, 0};
        } else {
            throw std::invalid_argument("Please choose either navierField or timeDependentNavierField if analyzing linear fields.");
        }
    }
}

/**
 * Return the indices of the contact point in 2D or the coordinates most closely matching the given parameter.
 * This function works differently for 2D and 3D. If the simulation is 2D, it ignores the input parameter and returns
 * the contact point by checking where the sign of the LevelSet field changes.  It uses a convex combination of
 * neighboring LevelSet values depending on whether the right or the left contact point is being tracked.
 *
 * In 3D, it returns the coordinates of the point according to the reference ODE. In other words, it is assumed that in 3D,
 * the solution of the ODE exactly matches the actual contact point.
 *
 * @param point A point in the space of the Level set field.
 * @param indexOnly Whether to return only the indices of the point or its actual coordinates. Default is false.
 * @return The contact point
 */
Vector LevelSet::getContactPoint(int timestep, bool indexOnly /* = false */) const {
    if (numZ == 1) {
        if (trackedCP == "left") {
            double initSign = this->at(0, 0, 0)/std::abs(this->at(0, 0, 0));

            if (initSign < 0)
                throw std::runtime_error("CP_Exit_Simulation_Plane");

            for (int i = 0; i < numX; i++) {
                if (this->at(i, 0, 0)*initSign < 0) {

                    if (indexOnly)
                        return {double(i) + 0.5, 0, 0}; //Add 0.5 to compensate floating-point errors
                    double alpha = this->at(i,0,0)-this->at(i - 1, 0, 0);
                    if(std::abs(alpha) < 1e-12) {
                        throw std::runtime_error("LevelSet::getContactPoint: \nDifference of LevelSet values at contact point too small for convex combination");
                    }

                    alpha = this->at(i,0,0)/(alpha);
                    return {(1 - alpha)*i*dx + alpha*(i - 1)*dx, 0, 0};
                }
            }
        } else {
            double initSign = this->at(numX - 1, 0, 0)/std::abs(this->at(numX - 1, 0, 0));

            if (initSign < 0)
                throw std::runtime_error("CP_Exit_Simulation_Plane");

            for (int i = numX - 1; i >= 0; i--) {
                if (this->at(i, 0, 0)*initSign < 0) {

                    if (indexOnly)
                        return {double(i) + 0.5, 0, 0};
                    double alpha = this->at(i + 1, 0, 0) - this->at(i, 0, 0);
                    if(std::abs(alpha) < 1e-12) {
                        throw std::runtime_error("LevelSet::getContactPoint: \nDifference of LevelSet values at contact point too small for convex combination");
                    }

                alpha = this->at(i + 1, 0, 0) / alpha;
                return {alpha*i*dx + (1 - alpha)*(i+1)*dx, 0, 0};
                }
            }
        }
        throw std::runtime_error("CP_Exit_Simulation_Plane");
    }

    else {
        const Vector ref = positionReference[timestep];
         if (ref[0] < 0 || ref[1] < 0 || ref[2] < 0 ||
             ref[0] > this->numX*dx || ref[1] > this->numY*dy || ref[2] > this->numZ*dz) {
                 throw std::runtime_error("CP_Exit_Simulation_Plane");
             }
        if (indexOnly) {
            Vector cell = {0, 0, 0};
            for (int i = 0; i < numX; ++i)
                for (int k = 0; k < numZ; ++k) {
                    Vector current = {i*dx, 0, k*dz};
                    Vector temp = {int(cell[0]) * dx, 0, int(cell[2]) * dz};
                    if ( abs(current - ref) < abs(temp - ref))
                        cell = {double(i) + 0.5, 0, double(k) + 0.5};
                }
            return cell;
        } else {
            return ref;
        }
    }
}

/**
 * Get the current indices of the contact point.
 *
 * This is a wrapper function for LevelSet::getContactPoint.
 * @param timestep The timestep
 * @return The indices of the contact point.
 */
array<int, 3> LevelSet::getContactPointIndices(int timestep) const {
    Vector vec = getContactPoint(timestep, true);
    return { int(vec[0]), int(vec[1]), int(vec[2]) };
}

/**
 * Calculate the normal vector.
 *
 * This function uses finite differences to calculate the normal vector of the Level set field at cell
 * @param cell The indices of the contact point.
 * @param cellIsCP Whether cell is the index-vector of the contact point (default is true).
 * @param normalizeVector Whether the resulting gradient vector is normalized (default is true).
 * @param findCPin2D Whether to take the passed coordinates and calculate the normal there or whether to find
 * the contact point automatically (using the change in sign). Only applicable in 2D.
 * @return The normal vector at cell.
 */
Vector LevelSet::getNormalVector(array<int, 3> cell, bool useInterpolation /* = true */, bool normalizeVector /* = true */, bool findCPin2D /* = true */) const {
    double normalX = 0, normalY = 0, normalZ = 0;

    // TODO Mathis
    //useInterpolation = false;

    if (numZ == 1 && findCPin2D && useInterpolation) {

        if (trackedCP == "left") {
             // find root of phi, alpha: coefficient for convex combination
             double alpha = this->at(cell[0],0,0)-this->at(cell[0]-1,0,0);

             if(std::abs(alpha) < 1E-12){
               throw std::runtime_error("LevelSet::getContactAngle: \nDifference of LevelSet values at contact point too small for convex combination");
             }

             alpha = this->at(cell[0],0,0)/(alpha);

            //Calculate angle at this cell with finite differences
            normalX = alpha*(this->at(cell[0], cell[1], cell[2]) - this->at(cell[0] - 2, cell[1], cell[2]))/(2*dx)
            + (1-alpha)*(this->at(cell[0] + 1, cell[1], cell[2]) - this->at(cell[0] - 1, cell[1], cell[2]))/(2*dx);
            normalY = alpha*(-this->at(cell[0] - 1, cell[1] + 2, cell[2])
                              + 4.0*this->at(cell[0] - 1, cell[1] + 1, cell[2])
                              - 3.0*this->at(cell[0] - 1, cell[1], cell[2]))/(2*dy)
                              + (1-alpha)*(-this->at(cell[0], cell[1] + 2, cell[2])
                              + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
                              - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dy); // second order difference quotient

        } else if (trackedCP == "right") {
            double alpha = this->at(cell[0] + 1, 0, 0) - this->at(cell[0], 0, 0);
            if(std::abs(alpha) < 1E-12) {
               throw std::runtime_error("LevelSet::getContactAngle: \nDifference of LevelSet values at contact point too small for convex combination");
            }

            alpha = this->at(cell[0] + 1, 0, 0) / alpha;

            normalX = alpha*(this->at(cell[0] + 1, 0, 0) - this->at(cell[0] - 1, 0, 0))/(2*dx)
                    + (1 - alpha)*(this->at(cell[0] + 2, 0, 0) - this->at(cell[0], 0, 0))/(2*dx);
            normalY = alpha*(-this->at(cell[0], cell[1] + 2, cell[2])
                              + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
                              - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dy)
                              + (1 - alpha)*(-this->at(cell[0] + 1, cell[1] + 2, cell[2])
                              + 4.0*this->at(cell[0] + 1, cell[1] + 1, cell[2])
                              - 3.0*this->at(cell[0] + 1, cell[1], cell[2]))/(2*dy);
        }

    } else {
        if (cell[0] > 0 && cell[0] < numX - 1) {
                normalX =    (this->at(cell[0] + 1, cell[1], cell[2])
                            - this->at(cell[0] - 1, cell[1], cell[2]))/(2*dx);
        } else if (cell[0] == 0) {
                normalX = (-1*this->at(cell[0] + 2, cell[1], cell[2])
                           +4*this->at(cell[0] + 1, cell[1], cell[2])
                           -3*this->at(cell[0], cell[1], cell[2]))/(2*dx);
        } else {
                normalX = (   this->at(cell[0] - 2, cell[1], cell[2])
                           -4*this->at(cell[0] - 1, cell[1], cell[2])
                           +3*this->at(cell[0], cell[1], cell[2]))/(2*dx);
        }

        if (cell[1] > 0 && cell[1] < numY - 1) {
            normalY =   ( this->at(cell[0], cell[1] + 1, cell[2])
                        - this->at(cell[0], cell[1] - 1, cell[2])) / (2 * dy);
        } else if (cell[1] == 0) {
            normalY = (-1*this->at(cell[0], cell[1] + 2, cell[2])
                       +4*this->at(cell[0], cell[1] + 1, cell[2])
                       -3*this->at(cell[0], cell[1], cell[2])) / (2 * dy);
        } else {
            normalY =    (this->at(cell[0], cell[1] - 2, cell[2])
                       -4*this->at(cell[0], cell[1] - 1, cell[2])
                       +3*this->at(cell[0], cell[1], cell[2])) / (2 * dy);
        }

        if (numZ > 1) {
            if (cell[2] > 0 && cell[2] < numZ - 1) {
                    normalZ =   (this->at(cell[0], cell[1], cell[2] + 1)
                               - this->at(cell[0], cell[1], cell[2] - 1))/(2*dz);
            } else if (cell[2] == 0) {
                    normalZ =(-1*this->at(cell[0], cell[1], cell[2] + 2)
                              +4*this->at(cell[0], cell[1], cell[2] + 1)
                              -3*this->at(cell[0], cell[1], cell[2]))/(2*dz);
            } else {
                    normalZ =   (this->at(cell[0], cell[1], cell[2] - 2)
                              -4*this->at(cell[0], cell[1], cell[2] - 1)
                              +3*this->at(cell[0], cell[1], cell[2]))/(2*dz);
            }
        }
    }

    Vector normal = {normalX, normalY, normalZ};

    //TODO Mathis: What if grad phi is numerically zero?
    if(abs(normal)<1E-15){
        throw std::runtime_error("Exit_Normal_Vector_0");
    }

    if (normalizeVector)
        normal = normal/abs(normal);
    return normal;
}

/**
 * @brief Get normal vector at given indices.
 *
 * This is a wrapper function for getNormalVector(array<int, 3>)
 * @param i Index along the x-axis
 * @param j Index along the y-axis
 * @param k Index along the z-axis
 * @return
 */
Vector LevelSet::getNormalVector(int i, int j, int k) const {
    array<int, 3> cell {i, j, k};
    return getNormalVector(cell);
}

/**
 * Calculate the tangential vector \f$\tau\f$.
 *
 * This function calculates \f$\tau\f$, where \f$\tau\f$ and the normal vector \f$\n_{\Sigma}\f$ define a plane perpendicular to the x-z plane.
 * @param normal The corresponding normal vector
 * @return The tangential vector
 */
Vector LevelSet::getTangentialVector(Vector normal) const {
//    double tau1, tau2, tau3;
    Vector tau;
    if ( std::abs(normal[1] - 1) > 1e-12) {
//        if (std::abs(normal[0]) < 1e-12) {
//            throw std::runtime_error("Error: Normal_x too small");
//        }
//        tau1 = 1;
//        tau2 = -(normal[2]*normal[2] / (normal[0]*normal[1]) + normal[0]/normal[1]);
//        tau3 = normal[2] / normal[0];   // TODO: What happens when normal[0] == 0 ? Check again.
        tau = cross(normal, {0, 0, 1});
    } else {
//        tau1 = 0;
//        tau2 = 1;
//        tau3 = 0;
        tau = cross(normal, {1, 0, 0});
    }

    // tau = {tau1, tau2, tau3};
    tau = tau/abs(tau);

    if (tau[1] > 0)
        return tau;
    else
        return -1*tau;
}

/**
 * Calculate the contact angle.
 * This function uses finite differences to calculate the normal vector of the Level set field at cell and
 * thus calculate the contact angle. This is a wrapper function for LevelSet::getNormalVector.
 *
 * @param timestep The timestep when the contact angle is calculated.
 * @return The contact angle in degrees
 */
double LevelSet::getContactAngleInterpolated(int timestep) {
    Vector normal;
    array<int, 3> cell = getContactPointIndices(timestep);
    if (numZ == 1)
        normal = getNormalVector(cell);
    else {
        Vector contactPoint = positionReference[timestep];
        double x = contactPoint[0];
        double y = contactPoint[1];
        double z = contactPoint[2];
        int i = floor(x/dx);
        int j = ceil(y/dy);
        int k = floor(z/dz);

        double lambda = x/dx - i;
        double mu = z/dz - k;

        normal = (1 - mu) * ((1 - lambda) * getNormalVector(i, j, k) + lambda*getNormalVector(i + 1, j, k))
           + mu * ((1 - lambda) * getNormalVector(i, j, k + 1) + lambda*getNormalVector(i + 1, j, k + 1));
    }

    return acos(normal[1]);
}

/**
 * Calculate the contact angle at a given time using the explicit Euler method.
 *
 * @param dt The length of a timestep
 * @param timestep The index of the timestep
 * @param n_sigma_init The initial normal vector of the interface
 * @param CP The current position of the contact point
 */
void LevelSet::referenceNormalExplicitEuler(double dt, int last_timestep, Vector normal_init) {
    Vector &n_sigma = normal_init;
	Vector deriv {0, 0, 0};
	for (int i = 0; i < last_timestep; i++) {
        normalReference[i] = n_sigma;
        angleReference[i] = acos(n_sigma[1])/M_PI*180;
	    Vector CP = positionReference[i];
		deriv = -1*transpose(field->gradAt(i*dt, CP[0], CP[1], CP[2]))*n_sigma
                + ((field->gradAt(i*dt, CP[0], CP[1], CP[2])*n_sigma)*n_sigma)*n_sigma;
		n_sigma = deriv*dt + n_sigma;
		n_sigma = n_sigma/abs(n_sigma);
	}
}

Vector LevelSet::referenceNormalExplicitEulerSingle(double dt, int last_timestep, std::vector<Vector> backwardsPoints) {
    Vector n_sigma = expectedNormalVector(backwardsPoints[backwardsPoints.size() - 1]);
	Vector deriv {0, 0, 0};
	for (int i = 0; i < last_timestep; i++) {
	    Vector CP = backwardsPoints[backwardsPoints.size() - i - 1];
		deriv = -1*transpose(field->gradAt(i*dt, CP[0], CP[1], CP[2]))*n_sigma
                + ((field->gradAt(i*dt, CP[0], CP[1], CP[2])*n_sigma)*n_sigma)*n_sigma;
		n_sigma = deriv*dt + n_sigma;
		n_sigma = n_sigma/abs(n_sigma);
	}
    return n_sigma;
}


/**
 * Calculate the contact angle at a given time for the navier field using the analytic solution.
 *
 * @param t The time at which to calculate the contact angle
 * @param c1 A parameter of the navier field
 * @param c2 A parameter of the navier field
 * @param theta0 The contact angle at time t == 0
 * @return The reference contact angle in degrees
 */
void LevelSet::referenceAngleLinearField(double dt, int last_timestep, double theta0) {
    double c1 = field->getC1();
    double c2 = field->getC2();
    for (int i = 0; i < last_timestep; ++i) {
        double t = i*dt;
        double temp;
        if (field->getName() == "navierField")
            if (trackedCP == "left") {
                temp = M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*t) - c2 * (exp(2*c1*t) - 1)/(2*c1));
            } else {
                temp =  M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*t) + c2 * (exp(2*c1*t) - 1)/(2*c1));
            }
        else if (field->getName() == "timeDependentNavierField") {
            double tau = field->getTau();
            if (trackedCP == "left") {
                temp = M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) - c2 * (exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) - 1)/(2*c1));
            } else {
                temp = M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) + c2 * (exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) - 1)/(2*c1));
            }
        }
        else
            throw std::invalid_argument("Please choose either navierField or timeDependentNavierField if analyzing linear fields.");

        angleReference[i] = temp/M_PI*180;
    }
}

/**
 * Calculate the reference curvature at the contact point with the explicit Euler method at a given time.
 * Currently, there are two different implementations in the 3D case: One based on the field-formulation of the normal vector,
 * and another (newer) based on sectional curvatures. TODO: Test both and decide which to keep.
 *
 * @param dt The length of a timestep
 * @param last_timestep The index of the last timestep.
 * @param initCurvature The initial curvature at time t == 0
 * @param CP The coordinates of the contact point
 * @param cell The indices of the contact point
 */
void
LevelSet::referenceCurvatureExplicitEuler(double dt, int last_timestep, double initCurvature, std::string reference3D) {

    if (numZ == 1) {
        double curvature = initCurvature;
        for (int i = 0; i < last_timestep; i++) {
            curvatureReference[i] = curvature;
            Vector CP = positionReference[i];
            Vector normal = normalReference[i];
            Vector tau = getTangentialVector(normal);
            // temp is the second derivative of v in the tau direction
            Vector temp = field->secondPartial(i*dt, CP[0], CP[1], CP[2], tau);
            curvature = curvature + dt*(temp*normal - 3*curvature*(field->gradAt(i*dt, CP[0], CP[1], CP[2])*tau)*tau);
        }
    } else if (reference3D == "fieldFormulation") {
        double d = std::min(std::min(dx, dy), dz);

        for (int i = 0; i < last_timestep; i++) {
            Vector CP = positionReference[i];
            int local = 2;
            int sidelengthZ;
            if (numZ > 1)
                sidelengthZ = 2*local + 1;
            else
                sidelengthZ = 1;

            // Declare a field of normal vectors
            Field<Vector> localField(2*local + 1, local + 1, sidelengthZ, dx, dy, dz);

            for (int x = -local; x <= local; x++)
                for (int y = 0; y <= local; y++)
                    for (int z = -sidelengthZ/2; z <= sidelengthZ/2; z++) {
                        Vector temp = CP + Vector({d*x, d*y, d*z});

                        if (temp[0] < 0 || temp[0] > numX*dx || temp[1] < 0 || temp[1] > numY*dy || temp[2] < 0 || temp[2] > numZ*dz)
                            continue;

                        std::vector<Vector> backwardsPoints = referenceContactPointBackwards(dt, i, temp);   // Backwards in time - Cha
                        Vector normal = referenceNormalExplicitEulerSingle(dt, i, backwardsPoints); // Forwards in time - Cha.  Cha-Cha!

                        if (this->numZ > 1)
                            localField.at(x+local, y, z+local) = normal;
                        else
                            localField.at(x+local, y, 0) = normal;
                    }

            // Calculate divergence of localField at cell, which by definition
            // is in the "middle" of localField
            double dnx_dx, dnx_dy, dnx_dz, dny_dx, dny_dy, dny_dz, dnz_dx, dnz_dy, dnz_dz;
            dnx_dx = dnx_dy = dnx_dz = dny_dx = dny_dy = dny_dz = dnz_dx = dnz_dy = dnz_dz = 0;


            //TODO : FIXME
            std::array<int, 3> cell = {numX/2, numY/2, numZ/2 };

            if (cell[0] == 0) {
                dnx_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[0]
                        + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[0]
                        - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*d);

                dny_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[1]
                        + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[1]
                        - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*d);

                if (numZ > 1) {
                    dnz_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[2]
                            + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[2]
                            - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*d);
                }

            } else if (cell[0] == numX - 1) {
                dnx_dx = (localField.at(local - 2, 0, sidelengthZ/2)[0]
                        - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[0]
                        + 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*d);

                dny_dx = (localField.at(local - 2, 0, sidelengthZ/2)[1]
                        - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[1]
                        + 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*d);

                if (numZ > 1) {
                    dnz_dx = (localField.at(local - 2, 0, sidelengthZ/2)[2]
                            - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[2]
                            + 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*d);
                }
            } else {
                dnx_dx = (localField.at(local + 1, 0, sidelengthZ/2)[0] - localField.at(local - 1, 0, sidelengthZ/2)[0]) / (2*d);
                dny_dx = (localField.at(local + 1, 0, sidelengthZ/2)[1] - localField.at(local - 1, 0, sidelengthZ/2)[1]) / (2*d);
                dnz_dx = (localField.at(local + 1, 0, sidelengthZ/2)[2] - localField.at(local - 1, 0, sidelengthZ/2)[2]) / (2*d);
            }

            dnx_dy = (-localField.at(local, 2, sidelengthZ/2)[0]
                    + 4.0*localField.at(local, 1, sidelengthZ/2)[0]
                    - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*d);


            dny_dy = (-localField.at(local, 2, sidelengthZ/2)[1]
                    + 4.0*localField.at(local, 1, sidelengthZ/2)[1]
                    - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*d);



            if (numZ > 1) {

                if (cell[2] == 0) {
                    dnx_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[0]
                            + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[0]
                            - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*d);
                    dny_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[1]
                            + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[1]
                            - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*d);
                    dnz_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[2]
                            + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[2]
                            - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*d);
                } else if (cell[2] == numZ - 1) {
                    dnx_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[0]
                            - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[0]
                            + 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*d);
                    dny_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[1]
                            - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[1]
                            + 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*d);
                    dnz_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[2]
                            - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[2]
                            + 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*d);
                } else {
                    dnx_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[0] - localField.at(local, 0, sidelengthZ/2 - 1)[0]) / (2*d);
                    dny_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[1] - localField.at(local, 0, sidelengthZ/2 - 1)[1]) / (2*d);
                    dnz_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[2] - localField.at(local, 0, sidelengthZ/2 - 1)[2]) / (2*d);
                }

                dnz_dy = (-localField.at(local, 2, sidelengthZ/2)[2]
                            + 4.0*localField.at(local, 1, sidelengthZ/2)[2]
                            - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*d);
            }

            Vector normal = localField.at(local, 0, sidelengthZ/2);
            Vector tau1 = getTangentialVector(normal);
            Vector tau2 = cross(normal, tau1);

            Vector row1 = {dnx_dx, dnx_dy, dnx_dz};
            Vector row2 = {dny_dx, dny_dy, dny_dz};
            Vector row3 = {dnz_dx, dnz_dy, dnz_dz};

            Matrix gradNormal = {row1, row2, row3};

            double kappa = -1*(gradNormal*tau1)*tau1 - 1*(gradNormal*tau2)*tau2;

            curvatureReference[i] = kappa;
        }
    } else if (reference3D == "sectionalCurvature") {
        double curvature = initCurvature;
        double secCurvA = 0.5*curvature;
        double secCurvB = 0.5*curvature;
        double secCurvC = 0.5*curvature;
        for (int i = 0; i < last_timestep; i++) {
            curvatureReference[i] = curvature;
            sectionalCurvatureAReference[i] = secCurvA;
            sectionalCurvatureBReference[i] = secCurvB;
            sectionalCurvatureCReference[i] = secCurvC;
            Vector CP = positionReference[i];
            Vector normal = normalReference[i];
            Vector tangentA = tangentAReference[i];
            Vector tangentB = tangentBReference[i];
            Vector tangentC = tangentCReference[i];
            // temp is the second derivative of v in the tau direction
            Vector tempA = field->secondPartial(i*dt, CP[0], CP[1], CP[2], tangentA);
//            curvature = curvature + dt*(temp*normal - 3*curvature*(field->gradAt(i*dt, CP[0], CP[1], CP[2])*tau)*tau);

            // TODO: Remove debug variables:
            secCurvA = secCurvA + dt * (tempA*normal + secCurvA * ((field->gradAt(i*dt, CP[0], CP[1], CP[2]))*normal)*normal
                    - 2*secCurvA*(field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentA)*tangentA);

            Vector tempB = field->secondPartial(i*dt, CP[0], CP[1], CP[2], tangentB);
            secCurvB = secCurvB + dt * (tempB*normal + secCurvB * ((field->gradAt(i*dt, CP[0], CP[1], CP[2]))*normal)*normal
                    - 2*secCurvB*(field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentB)*tangentB);

            Vector tempC = field->secondPartial(i*dt, CP[0], CP[1], CP[2], tangentC);
            secCurvC = secCurvC + dt * (tempC*normal + secCurvC * ((field->gradAt(i*dt, CP[0], CP[1], CP[2]))*normal)*normal
                    - 2*secCurvC*(field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentC)*tangentC);

            double mixedElement = (secCurvC - alpha1[i]*alpha1[i]*secCurvA - alpha2[i]*alpha2[i]*secCurvB) / (2*alpha1[i]*alpha2[i]);
            double projection = tangentA*tangentB;
            Eigen::Matrix2f tau1mat;
            tau1mat << 1, projection, projection, 1;   // TODO: Einlesen von Werten checken
            Eigen::Vector2f rightHandTau1;
            rightHandTau1 << 1, 0;
            Eigen::Vector2f c1andc2 = tau1mat.colPivHouseholderQr().solve(rightHandTau1);
            double c1 = c1andc2[0];
            double c2 = c1andc2[1];
            Eigen::Matrix2f tau2mat;
            tau2mat << 1, projection, projection, 1;
            Eigen::Vector2f rightHandTau2;
            rightHandTau2 << 0, 1;
            Eigen::Vector2f c3andc4 = tau2mat.colPivHouseholderQr().solve(rightHandTau2);
            double c3 = c3andc4[0];
            double c4 = c3andc4[1];

            curvature = c1*secCurvA + c4*secCurvB + (c2+c3)*mixedElement;
        }
    } else {
        throw std::invalid_argument("No valid curvature 3D reference selected.\n");
    }
}

void LevelSet::referenceTangentExplicitEuler(double dt, int last_timestep, Vector normal_init) {
    Vector tangentA = getTangentialVector(normal_init);
    Vector tangentB = cross(normal_init, tangentA);
    Vector tangentC = tangentA + tangentB;
    tangentC = tangentC / abs(tangentC);

    Vector derivA {0, 0, 0};
    Vector derivB {0, 0, 0};
    Vector derivC {0, 0, 0};
    for (int i = 0; i < last_timestep; i++) {
        tangentAReference[i] = tangentA;
        tangentBReference[i] = tangentB;
        tangentCReference[i] = tangentC;
        Eigen::Matrix<float, 3, 2> Mat;
        Mat << tangentA[0], tangentB[0], tangentA[1], tangentB[1], tangentA[2], tangentB[2];
        Eigen::Vector3f tangentCEigen;
        tangentCEigen << tangentC[0], tangentC[1], tangentC[2];
        Eigen::Vector2f alphas = Mat.colPivHouseholderQr().solve(tangentCEigen);
        alpha1[i] = alphas[0];
        alpha2[i] = alphas[1];
        Vector CP = positionReference[i];
        // This is the formula copied from the normal reference. TODO: Delete later.
//        deriv = -1 * transpose(field->gradAt(i*dt, CP[0], CP[1], CP[2])) * tangent
//                + ((field->gradAt(i*dt, CP[0], CP[1], CP[2]) * tangent) * tangent) * tangent;
        derivA = field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentA
                - ((field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentA)*tangentA)*tangentA;
        derivB = field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentB
                - ((field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentB)*tangentB)*tangentB;
        derivC = field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentC
                - ((field->gradAt(i*dt, CP[0], CP[1], CP[2])*tangentC)*tangentC)*tangentC;

        tangentA = derivA * dt + tangentA;
        tangentA = tangentA / abs(tangentA);
        tangentB = derivB * dt + tangentB;
        tangentB = tangentB / abs(tangentB);
        tangentC = derivC * dt + tangentC;
        tangentC = tangentC / abs(tangentC);
    }
}


/**
 * @brief LevelSet::referenceCurvatureDeriv3D Calculate the curvature derivative at time t = 0 in 3D.
 * @param initCurvature The initial curvature. (= -1/R)
 */
double LevelSet::referenceCurvatureDeriv3D(double initCurvature, Matrix expNormalVectorGrad) {
    Vector CP = positionReference[0];
    Vector normal = normalReference[0];
    Vector tau0 = getTangentialVector(normal);
    Vector tau1 = cross(tau0, normal);
    std::array<Vector, 2> tau = {tau0, tau1};

    double ret = 0;
    for (int i = 0; i < 2; i++) {
        ret += field->secondPartial(0, CP[0], CP[1], CP[2], tau[i])*normal;
        for (int j = 0; j < 2; j++)
            ret += 2 * (tau[j] * (field->gradAt(0, CP[0], CP[1], CP[2])* tau[i])) * (tau[j] * (expNormalVectorGrad * tau[i]));
    }

    ret += initCurvature * (transpose(field->gradAt(0, CP[0], CP[1], CP[2]))*normal) * normal;

    return ret;
}

/**
 * Calculate the reference curvature for a specific case of the navier field (c0 = 0 = c2)
 *
 * @param t The time
 * @param init_curvature The initial curvature
 */
void LevelSet::referenceCurvatureLinearField(double dt, int timesteps, double initCurvature) {
    double curvature = initCurvature;
    double c1 = field->getC1();
    for (int i = 0; i < timesteps; ++i) {
        double t = dt*i;
        curvature = initCurvature * exp(3*c1*t);
        curvatureReference[i] = curvature;
    }
}

/**
 * Calculate the reference curvature for a specific case of the quadratic (c2 = 0)
 *
 * @param t The time
 * @param init_curvature The initial curvature
 */
void LevelSet::referenceCurvatureQuadraticField(double dt, int timesteps, double initCurvature) {
    double curvature = initCurvature;
    for (int i = 0; i < timesteps; ++i) {
        double c1 = field->getC1();
        double c3 = field->getC3();
        double t = dt*i;
        curvature = initCurvature * exp(3*c1* t) + (2.0/3) * (c3/c1) * (1 - exp(3 *c1* t));
        curvatureReference[i] = curvature;
    }
}

/**
 * Calculate the curvature of the droplet at the contact point.
 *
 * First, it calculates a local field of normal vectors around the contact point with dimension 5x3x5 cells in 2D
 * and 5x3x1 in 2D.
 * Afterwards, it uses the local field to calculate the divergence in the middle of the slice with y = 0 and uses this
 * to calculate the curvature.
 *
 * @param cell The indices of the contact point
 * @return The curvature
 */
double LevelSet::getCurvature(array<int, 3> cell) const {
    /* Define what number of cells in each direction(excluding the main cell) are considered local.
       The resulting normal vector field is defined on a cuboid with sidelength 2*local + 1 */
    int local = 2;
    int sidelengthZ;
    if (numZ > 1)
        sidelengthZ = 2*local + 1;
    else
        sidelengthZ = 1;

    // Declare a field of normal vectors
    Field<Vector> localField(2*local + 1, local + 1, sidelengthZ, dx, dy, dz);

    for (int x = -local; x <= local; x++)
        for (int y = 0; y <= local; y++)
            for (int z = -sidelengthZ/2; z <= sidelengthZ/2; z++) {
                array<int, 3> temp = {x, y, z};
                temp = temp + cell;

                if (temp[0] < 0 || temp[0] > numX - 1 || temp[1] < 0 || temp[1] > numY - 1 || temp[2] < 0 || temp[2] > numZ - 1)
                    continue;

                Vector normal = getNormalVector(temp, false);

                if (this->numZ > 1)
                    localField.at(x+local, y, z+local) = normal;
                else
                    localField.at(x+local, y, 0) = normal;
            }

    // Calculate divergence of localField at cell, which by definition
    // is in the "middle" of localField
    double dnx_dx, dnx_dy, dnx_dz, dny_dx, dny_dy, dny_dz, dnz_dx, dnz_dy, dnz_dz;
    dnx_dx = dnx_dy = dnx_dz = dny_dx = dny_dy = dny_dz = dnz_dx = dnz_dy = dnz_dz = 0;

    if (cell[0] == 0) {
        dnx_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[0]
                + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[0]
                - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dx);

        dny_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[1]
                + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[1]
                - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dx);

        if (numZ > 1) {
            dnz_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[2]
                    + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[2]
                    - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dx);
        }

    } else if (cell[0] == numX - 1) {
        dnx_dx = (localField.at(local - 2, 0, sidelengthZ/2)[0]
                - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[0]
                + 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dx);

        dny_dx = (localField.at(local - 2, 0, sidelengthZ/2)[1]
                - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[1]
                + 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dx);

        if (numZ > 1) {
            dnz_dx = (localField.at(local - 2, 0, sidelengthZ/2)[2]
                    - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[2]
                    + 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dx);
        }
    } else {
        dnx_dx = (localField.at(local + 1, 0, sidelengthZ/2)[0] - localField.at(local - 1, 0, sidelengthZ/2)[0]) / (2*dx);
        dny_dx = (localField.at(local + 1, 0, sidelengthZ/2)[1] - localField.at(local - 1, 0, sidelengthZ/2)[1]) / (2*dx);
        dnz_dx = (localField.at(local + 1, 0, sidelengthZ/2)[2] - localField.at(local - 1, 0, sidelengthZ/2)[2]) / (2*dx);
    }

    dnx_dy = (-localField.at(local, 2, sidelengthZ/2)[0]
            + 4.0*localField.at(local, 1, sidelengthZ/2)[0]
            - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dy);


    dny_dy = (-localField.at(local, 2, sidelengthZ/2)[1]
            + 4.0*localField.at(local, 1, sidelengthZ/2)[1]
            - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dy);



    if (numZ > 1) {

        if (cell[2] == 0) {
            dnx_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[0]
                     + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[0]
                     - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dz);
            dny_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[1]
                     + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[1]
                     - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dz);
            dnz_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[2]
                     + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[2]
                     - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dz);
        } else if (cell[2] == numZ - 1) {
            dnx_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[0]
                    - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[0]
                    + 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dz);
            dny_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[1]
                    - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[1]
                    + 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dz);
            dnz_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[2]
                    - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[2]
                    + 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dz);
        } else {
            dnx_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[0] - localField.at(local, 0, sidelengthZ/2 - 1)[0]) / (2*dz);
            dny_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[1] - localField.at(local, 0, sidelengthZ/2 - 1)[1]) / (2*dz);
            dnz_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[2] - localField.at(local, 0, sidelengthZ/2 - 1)[2]) / (2*dz);
        }

        dnz_dy = (-localField.at(local, 2, sidelengthZ/2)[2]
                    + 4.0*localField.at(local, 1, sidelengthZ/2)[2]
                    - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dy);
    }

    Vector normal = localField.at(local, 0, sidelengthZ/2);
    Vector tau1 = getTangentialVector(normal);
    Vector tau2 = cross(normal, tau1);

    Vector row1 = {dnx_dx, dnx_dy, dnx_dz};
    Vector row2 = {dny_dx, dny_dy, dny_dz};
    Vector row3 = {dnz_dx, dnz_dy, dnz_dz};

    Matrix gradNormal = {row1, row2, row3};

    double kappa = -1*(gradNormal*tau1)*tau1 - 1*(gradNormal*tau2)*tau2;

    return kappa;
}

double LevelSet::getSectionalCurvature(array<int, 3> cell, Vector tau) {
// TODO: This function is (with exception of the last few lines) an exact copy-paste of getCurvature. TODO fix.

    /* Define what number of cells in each direction(excluding the main cell) are considered local.
       The resulting normal vector field is defined on a cuboid with sidelength 2*local + 1 */
    int local = 2;
    int sidelengthZ;
    if (numZ > 1)
        sidelengthZ = 2*local + 1;
    else
        sidelengthZ = 1;

    // Declare a field of normal vectors
    Field<Vector> localField(2*local + 1, local + 1, sidelengthZ, dx, dy, dz);

    for (int x = -local; x <= local; x++)
        for (int y = 0; y <= local; y++)
            for (int z = -sidelengthZ/2; z <= sidelengthZ/2; z++) {
                array<int, 3> temp = {x, y, z};
                temp = temp + cell;

                if (temp[0] < 0 || temp[0] > numX - 1 || temp[1] < 0 || temp[1] > numY - 1 || temp[2] < 0 || temp[2] > numZ - 1)
                    continue;

                Vector normal = getNormalVector(temp, false);

                if (this->numZ > 1)
                    localField.at(x+local, y, z+local) = normal;
                else
                    localField.at(x+local, y, 0) = normal;
            }

    // Calculate divergence of localField at cell, which by definition
    // is in the "middle" of localField
    double dnx_dx, dnx_dy, dnx_dz, dny_dx, dny_dy, dny_dz, dnz_dx, dnz_dy, dnz_dz;
    dnx_dx = dnx_dy = dnx_dz = dny_dx = dny_dy = dny_dz = dnz_dx = dnz_dy = dnz_dz = 0;

    if (cell[0] == 0) {
        dnx_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[0]
                  + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[0]
                  - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dx);

        dny_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[1]
                  + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[1]
                  - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dx);

        if (numZ > 1) {
            dnz_dx = (-localField.at(local + 2, 0, sidelengthZ/2)[2]
                      + 4.0*localField.at(local + 1, 0, sidelengthZ/2)[2]
                      - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dx);
        }

    } else if (cell[0] == numX - 1) {
        dnx_dx = (localField.at(local - 2, 0, sidelengthZ/2)[0]
                  - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[0]
                  + 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dx);

        dny_dx = (localField.at(local - 2, 0, sidelengthZ/2)[1]
                  - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[1]
                  + 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dx);

        if (numZ > 1) {
            dnz_dx = (localField.at(local - 2, 0, sidelengthZ/2)[2]
                      - 4.0*localField.at(local - 1, 0, sidelengthZ/2)[2]
                      + 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dx);
        }
    } else {
        dnx_dx = (localField.at(local + 1, 0, sidelengthZ/2)[0] - localField.at(local - 1, 0, sidelengthZ/2)[0]) / (2*dx);
        dny_dx = (localField.at(local + 1, 0, sidelengthZ/2)[1] - localField.at(local - 1, 0, sidelengthZ/2)[1]) / (2*dx);
        dnz_dx = (localField.at(local + 1, 0, sidelengthZ/2)[2] - localField.at(local - 1, 0, sidelengthZ/2)[2]) / (2*dx);
    }

    dnx_dy = (-localField.at(local, 2, sidelengthZ/2)[0]
              + 4.0*localField.at(local, 1, sidelengthZ/2)[0]
              - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dy);


    dny_dy = (-localField.at(local, 2, sidelengthZ/2)[1]
              + 4.0*localField.at(local, 1, sidelengthZ/2)[1]
              - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dy);



    if (numZ > 1) {

        if (cell[2] == 0) {
            dnx_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[0]
                      + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[0]
                      - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dz);
            dny_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[1]
                      + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[1]
                      - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dz);
            dnz_dz = (-localField.at(local, 0, sidelengthZ/2 + 2)[2]
                      + 4.0*localField.at(local, 0, sidelengthZ/2 + 1)[2]
                      - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dz);
        } else if (cell[2] == numZ - 1) {
            dnx_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[0]
                      - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[0]
                      + 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dz);
            dny_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[1]
                      - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[1]
                      + 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dz);
            dnz_dz = (localField.at(local, 0, sidelengthZ/2 - 2)[2]
                      - 4.0*localField.at(local, 0, sidelengthZ/2 - 1)[2]
                      + 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dz);
        } else {
            dnx_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[0] - localField.at(local, 0, sidelengthZ/2 - 1)[0]) / (2*dz);
            dny_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[1] - localField.at(local, 0, sidelengthZ/2 - 1)[1]) / (2*dz);
            dnz_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[2] - localField.at(local, 0, sidelengthZ/2 - 1)[2]) / (2*dz);
        }

        dnz_dy = (-localField.at(local, 2, sidelengthZ/2)[2]
                  + 4.0*localField.at(local, 1, sidelengthZ/2)[2]
                  - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dy);
    }
    Vector row1 = {dnx_dx, dnx_dy, dnx_dz};
    Vector row2 = {dny_dx, dny_dy, dny_dz};
    Vector row3 = {dnz_dx, dnz_dy, dnz_dz};

    Matrix gradNormal = {row1, row2, row3};

    double testVar = -1*(gradNormal*tau)*tau;

    return -1*(gradNormal*tau)*tau;
}

double LevelSet::getSectionalCurvatureInterpolated(int timestep, Vector tau) {
        Vector contactPoint = positionReference[timestep];
        double x = contactPoint[0];
        double z = contactPoint[2];
        int i = floor(x/dx);
        int j = 0;
        int k = floor(z/dz);

        double lambda = x/dx - i;
        double mu = z/dz - k;

//        Removed interpolation, as it seems to be causing problems.
        double curvature = (1 - mu) * ((1 - lambda) * getSectionalCurvature({i, j, k}, tau) + lambda*getSectionalCurvature({i + 1, j, k}, tau))
                    + mu * ((1 - lambda) * getSectionalCurvature({i, j, k + 1}, tau) + lambda* getSectionalCurvature({i + 1, j, k + 1}, tau));
//        double curvature = getSectionalCurvature({i, j, k}, tau);

        return curvature;
}

double LevelSet::getCurvature(int i, int j, int k) const {
    array<int, 3> cell{i, j, k};
    return getCurvature(cell);
}

double LevelSet::getCurvatureInterpolated(int timestep) const {

    double curvature = 0;

    if (numZ == 1) {
        array<int, 3> cell = getContactPointIndices(timestep);
        if (trackedCP == "left") {
             // find root of phi, alpha: coefficient for convex combination
             double alpha = this->at(cell[0],0,0)-this->at(cell[0]-1,0,0);

             if(std::abs(alpha) < 1E-12){
               throw std::runtime_error("LevelSet::getCurvatureInterpolated: \nDifference of LevelSet values at contact point too small for convex combination");
             }

             alpha = this->at(cell[0],0,0)/(alpha);

            //Calculate angle at this cell with finite differences
            /*normalX = alpha*(this->at(cell[0], cell[1], cell[2]) - this->at(cell[0]-2, cell[1], cell[2]))/(2*dx)
            + (1-alpha)*(this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2]))/(2*dx);
            normalY = alpha*(-this->at(cell[0]-1, cell[1] + 2, cell[2])
                              + 4.0*this->at(cell[0]-1, cell[1] + 1, cell[2])
                              - 3.0*this->at(cell[0]-1, cell[1], cell[2]))/(2*dy)
                              + (1-alpha)*(-this->at(cell[0], cell[1] + 2, cell[2])
                              + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
                              - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dy);*/ // second order difference quotient

             curvature = alpha * getCurvature({cell[0] - 1, cell[1], cell[2]}) + (1 - alpha) * getCurvature(cell);



        } else if (trackedCP == "right") {
            double alpha = this->at(cell[0] + 1, 0, 0) - this->at(cell[0], 0, 0);
            if(std::abs(alpha) < 1E-12) {
               throw std::runtime_error("LevelSet::getCurvatureInterpolated: \nDifference of LevelSet values at contact point too small for convex combination");
            }

            alpha = this->at(cell[0] + 1, 0, 0) / alpha;

//            normalX = alpha*(this->at(cell[0] + 1, 0, 0) - this->at(cell[0] - 1, 0, 0))/(2*dx)
//                    + (1 - alpha)*(this->at(cell[0] + 2, 0, 0) - this->at(cell[0], 0, 0))/(2*dx);
//            normalY = alpha*(-this->at(cell[0], cell[1] + 2, cell[2])
//                              + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
//                              - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dy)
//                              + (1 - alpha)*(-this->at(cell[0] + 1, cell[1] + 2, cell[2])
//                              + 4.0*this->at(cell[0] + 1, cell[1] + 1, cell[2])
//                              - 3.0*this->at(cell[0] + 1, cell[1], cell[2]))/(2*dy);

            curvature = alpha * getCurvature(cell) + (1 - alpha) * getCurvature({cell[0] + 1, cell[1], cell[2]});
        }
    } else {
        Vector contactPoint = positionReference[timestep];
        double x = contactPoint[0];
        double z = contactPoint[2];
        int i = floor(x/dx);
        int j = 0;
        int k = floor(z/dz);

        double lambda = x/dx - i;
        double mu = z/dz - k;

        curvature = (1 - mu) * ((1 - lambda) * getCurvature(i, j, k) + lambda*getCurvature(i + 1, j, k))
           + mu * ((1 - lambda) * getCurvature(i, j, k + 1) + lambda*getCurvature(i + 1, j, k + 1));
    }

    return curvature;
}

/**
 * Write the files necessary to visualize the field and its velocity field in Paraview.
 *
 * Write the XMF file to disk, as well as the grid of the Level set field. This is done only once.
 * At each point in time, write the Level set field to disk and call VelocityField::writeToFile to write its values as well.
 *
 * @param dt The width of a timestep
 * @param timestep The index of the timestep
 * @param total_timesteps The total number of timesteps
 * @param total_writesteps The total number of steps that will be written to disk
 * @param mainXmfFile A pointer to the XMF file.
 */
void
LevelSet::writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *mainXmfFile,
                      std::ofstream *tauXmfFile) {
    int Npoints = numX*numY*numZ;
    double *pointCoordinates = new double[Npoints*3];
    double *pointPhiValues = new double[Npoints];
    double *pointRValues = new double[Npoints];

    int index = 0;

    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++) {
            for (int i = 0; i < numX; i++) {
                if (timestep == 0) {
                    pointCoordinates[index] = i*dx;
                    pointCoordinates[index +1] = j*dy;
                    pointCoordinates[index +2] = k*dz;
                    index += 3;
                }
                pointPhiValues[i + j*numX + k*numX*numY] = this->at(i, j, k);
                Vector normalVector;
                try {
                    normalVector = getNormalVector({i, j, k}, false, true, false);
                } catch (std::runtime_error& e) {
                    std::string str = e.what();
                    if (str.compare("Exit_Normal_Vector_0")) {
                        normalVector = {0, 0, 0};
                    }
                }
                pointRValues[i + j*numX + k*numX*numY] = (field->gradAt(dt*timestep, i*dx, j*dy, k*dz) * normalVector) * normalVector;
            }
        }

    // If it is the first iteration, create coordinate file
    if (timestep == 0) {
        *mainXmfFile   << "<?xml version=\"1.0\" ?>\n"
                       "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
                       "<!ENTITY Npoints \"" + std::to_string(Npoints) + "\">\n"
                       "<!ENTITY numX \"" + std::to_string(numX) + "\">\n"
                       "<!ENTITY numY \"" + std::to_string(numY) + "\">\n"
                       "<!ENTITY numZ \"" + std::to_string(numZ) + "\">\n";
        *tauXmfFile << "<?xml version=\"1.0\" ?>\n"
                       "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
                       "<!ENTITY Npoints \"1\">\n"
                       "<!ENTITY numX \"1\">\n"
                       "<!ENTITY numY \"1\">\n"
                       "<!ENTITY numZ \"1\">\n";

        std::string temp = "]>\n"
                         "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
                         "<Domain>\n"
                         "<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
                         "<Time TimeType=\"List\">\n"
                         "<DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\""+ std::to_string(total_writesteps) + "\">\n";

        *mainXmfFile << temp;
        *tauXmfFile << temp;

        for (int i = 0; i < total_writesteps; i++) {
            *mainXmfFile << std::to_string(i*((double)total_timesteps/total_writesteps)*dt) << " ";
            *tauXmfFile << std::to_string(i*((double)total_timesteps/total_writesteps)*dt) << " ";
        }
        *mainXmfFile <<"</DataItem>\n" << "</Time>\n";
        *tauXmfFile <<"</DataItem>\n" << "</Time>\n";

        //Write field coordinates into binary file
        FILE *fieldFile;
        temp = outputDirectory + "data/field.bin";
        fieldFile = fopen(temp.data(), "wb");
        fwrite(pointCoordinates, sizeof(double), Npoints*3, fieldFile);
        fclose(fieldFile);
    }

    FILE *PhiFile;
    std::string filename = outputDirectory + "data/Phi_t="+ std::to_string(timestep*dt)+".bin";
    PhiFile = fopen(filename.data(), "wb");
	fwrite(pointPhiValues, sizeof(double), Npoints, PhiFile);
	fclose(PhiFile);

    FILE *rFile;
    std::string filenameRFile = outputDirectory + "data/r_t="+ std::to_string(timestep*dt)+".bin";
    rFile = fopen(filenameRFile.data(), "wb");
    fwrite(pointRValues, sizeof(double), Npoints, rFile);
    fclose(rFile);

    // TODO: Testen ob wirklich Reihenfolge numX,numY,numZ in NumberOfElements und Dimensions bei Topology und DataItem vertauscht ist

     *mainXmfFile << "<Grid GridType=\"Uniform\">\n"
             << "<Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"&numZ; &numY; &numX;\"/>\n"
             << "<Geometry GeometryType=\"ORIGIN_DXDYDZ\"> \n"
             << "<DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\">"
             << "0.00 0.00 0.00"
             << "</DataItem>"
             << "<DataItem Name=\"points\" NumberType=\"Float\" Dimensions=\"3\">\n"
             << std::to_string(dz) + " " + std::to_string(dy) + " " + std::to_string(dx)
             << "</DataItem>\n"
             << "</Geometry>\n"
             << "<Attribute Name =\"lvlset\" AttributeType=\"Scalar\" Center=\"Node\">\n"
             << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\"  Endian=\"Little\" Dimensions=\"&numZ; &numY; &numX;\">\n"
			 << "Phi_t=" + std::to_string(timestep*dt) +".bin\n"
			 << "</DataItem></Attribute>\n"
             << "<Attribute Name =\"sourceTermField\" AttributeType=\"Scalar\" Center=\"Node\">\n"
             << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\"  Endian=\"Little\" Dimensions=\"&numZ; &numY; &numX;\">\n"
             << "r_t=" + std::to_string(timestep*dt) +".bin\n"
             << "</DataItem></Attribute>\n"
             << "<Attribute Name =\"VelocityField\" AttributeType=\"Vector\" Center=\"Node\">\n"
             << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&numZ; &numY; &numX; 3\">\n";


    if (field->getName() == "timeDependentNavierField") {
		*mainXmfFile << "Vel_t=" + std::to_string(timestep * dt) + ".bin\n";
    }
    else {
		*mainXmfFile << "Vel_t=" + std::to_string(0 * dt) + ".bin\n";
    }

        *mainXmfFile << "</DataItem></Attribute>\n"
                 << "</Grid>\n";

    *tauXmfFile  << "<Grid>\n"
                 << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"1\"/>\n"
                 << "<Geometry GeometryType=\"XYZ\"> \n"
                 << "<DataItem Name=\"points\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
                 << "CP_reference_t=" + std::to_string(timestep*dt) +".bin\n"
                 << "</DataItem>\n</Geometry>\n"
                 << "<Attribute Name =\"Tangential Vector\" AttributeType=\"Vector\" Center=\"Cell\">\n"
                 << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
                 << "Tau_t=" + std::to_string(timestep*dt) +".bin\n"
                 << "</DataItem></Attribute>\n</Grid>\n";


    delete[] pointCoordinates;
    delete[] pointPhiValues;
    delete[] pointRValues;

    //Write velocity field
    if (field->getName() == "timeDependentNavierField")
    	field->writeToFile(timestep*dt);
    else if (timestep == 0)
    	field->writeToFile(0*dt);

    //Write tangential vector to file
    writeTangentialVectorToFile(dt, timestep);
}

/**
 * Writes the field of the tangential vector at the contact point at a given time to disk
 *
 * Writes the tangental vector field for visualization in Paraview. Since no XMF file is written,
 * simply calling this function is not enough for visualization. Thus, this function is called within
 * LevelSet::writeToFile, which does write an XMF file.
 *
 * @param t The time
 */
void LevelSet::writeTangentialVectorToFile(double dt, int timestep) {
    double *fieldTau = positionReference[timestep].data();
    Vector normal = getNormalVector(getContactPointIndices(timestep));
    Vector tangential = getTangentialVector(normal);

    std::string filenameTau = outputDirectory + "data/Tau_t=" + std::to_string(dt*timestep) + ".bin";
    std::string filenameFieldTau = outputDirectory + "data/CP_reference_t=" + std::to_string(dt*timestep) +".bin";
    FILE *tauFile;
    FILE *fieldTauFile;
    tauFile = fopen(filenameTau.data(), "wb");
    fieldTauFile = fopen(filenameFieldTau.data(), "wb");
    fwrite(tangential.data(), sizeof(double), 3, tauFile);
    fwrite(fieldTau, sizeof(double), 3, fieldTauFile);
    fclose(tauFile);
    fclose(fieldTauFile);
}

/**
 * Calculate the sum of the Level set field.
 * @return The value of the sum
 */
double LevelSet::sumLevelSet() {
    double temp = 0;
    for (int x = 0; x < numX; x++)
        for (int y = 0; y < numY; y++)
            for (int z = 0; z < numZ; z++)
                temp = temp + this->at(x, y, z);

    return temp;
}

/**
 * Initialize a droplet.
 * The droplet is initialized on the Level set field by setting the space occupied by the fluid to be negative,
 * the space occupied by the gas to be positive and the zero set to be the interface. More specifically, each
 * point on the grid is set to the value
 * \f[ \Phi = (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 \f]
 * Where \f$ x_0, y_0, z_0 \f$ are the coordinates of the droplet center and R is its radius.
 *
 * @param center The center of the droplet
 * @param radius The radius of the droplet
 */
void LevelSet::initSphere(Vector center, double radius) {

    if(numZ>1){

    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++)
            for (int i = 0; i < numX; i++)
                at(i, j, k) = 0.5*(pow(i*dx - center[0], 2)/radius + pow(j*dy - center[1], 2)/radius + pow(k*dz - center[2], 2)/radius - radius);
    }
    else if(numZ==1){
          for (int j = 0; j < numY; j++)
              for (int i = 0; i < numX; i++)
                  at(i, j, 0) = 0.5*(pow(i*dx - center[0], 2)/radius + pow(j*dy - center[1], 2)/radius - radius);

    }

}

/** Initialize a plane.
 *
 *  @param refPoint The reference point of the plane in degrees.
 *  @param angleA The angle between the initialized plane and the x-z-plane (= polar angle of normal vector) in deg.
 *  @param angleB The angle of rotation around the y-axis (= azimuthal angle) in deg.
 */
void LevelSet::initPlane(Vector refPoint, double polarAngle, double azimuthalAngle) {
    polarAngle = polarAngle/180*M_PI;
    azimuthalAngle = azimuthalAngle/180*M_PI;
    Vector normal = {sin(polarAngle) * cos(azimuthalAngle), cos(polarAngle), sin(polarAngle)* sin(azimuthalAngle)};

    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++)
            for (int i = 0; i < numX; i++)
                at(i, j, k) = normal * (Vector({i*dx, j*dy, k*dz}) - refPoint);

}

/** Initialize a paraboloid.
 * @param refPoint
 * @param stretchX Stretching factor in x-direction.
 * @param stretchZ Stretching factor in z-direction
 * @param heightMaximum Height of the maximum above the xz-plane.
 */
void LevelSet::initParaboloid(Vector refPoint, double stretchX, double stretchZ, double heightMaximum) {
    double x_0 = refPoint[0];
    double z_0 = refPoint[2];

    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++)
            for (int i = 0; i < numX; i++)
                at(i, j, k) = 0.5*stretchX * pow(i*dx - x_0, 2)  + 0.5*stretchZ* pow(k*dz - z_0, 2) + j*dy - heightMaximum;
}

/** Initialize an ellipsoid.
 * @param refPoint The reference point
 * @
 */
void LevelSet::initEllipsoid(Vector refPoint, double stretchX, double stretchY, double stretchZ) {
    double x_0 = refPoint[0];
    double y_0 = refPoint[1];
    double z_0 = refPoint[2];

    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++)
            for (int i = 0; i < numX; i++)
                at(i, j, k) = 0.5*stretchX * pow(i*dx - x_0, 2) + 0.5*stretchY * pow(j*dy - y_0, 2) + 0.5*stretchZ * pow(k*dz - z_0, 2) - 1;
}

Vector LevelSet::expectedNormalVector(Vector contactPoint) {
     std::vector<double>& params = shapeParams;
     Vector& refPoint = initCenter;

    Vector normal = {0, 0, 0};

    if (shape == InitShape::sphere) {
        normal = { contactPoint[0] - refPoint[0], contactPoint[1] - refPoint[1], contactPoint[2] - refPoint[2] };
    } else if (shape == InitShape::plane) {
        double polarAngle = params.at(0);
        double azimuthalAngle = params.at(1);
        normal =  { sin(polarAngle) * cos(azimuthalAngle), cos(polarAngle), sin(polarAngle)* sin(azimuthalAngle) };
    } else if (shape == InitShape::paraboloid) {
        double stretchX = params.at(0);
        double stretchZ = params.at(1);
        normal = { stretchX * (contactPoint[0] - refPoint[0]), 1, stretchZ * (contactPoint[2] - refPoint[2]) };
    } else if (shape == InitShape::ellipsoid) {
        double stretchX = params.at(0);
        double stretchY = params.at(1);
        double stretchZ = params.at(2);
        normal = { stretchX * (contactPoint[0] - refPoint[0]), stretchY * (contactPoint[1] - refPoint[1]), stretchZ * (contactPoint[2] - refPoint[2]) };
    }

    return normal/abs(normal);
}

Matrix LevelSet::expectedNormalVectorGradient(Vector contactPoint) {
    std::vector<double>& params = shapeParams;
    Vector vec0;
    Vector vec1;
    Vector vec2;

    double x = contactPoint[0];
    double y = 0;
    double z = contactPoint[2];
    double x0 = initCenter[0];
    double y0 = initCenter[1];
    double z0 = initCenter[2];


    if (shape == InitShape::sphere) {
        double n_abs = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
        double n_abs_3 = pow(n_abs, 3);
        vec0 = { (pow(y - y0, 2) + pow(z - z0, 2)) / n_abs_3, -(x - x0) * (y - y0) / n_abs_3, - (x - x0) * (z - z0) / n_abs_3 };
        vec1 = { -(x - x0) * (y - y0) / n_abs_3, (pow(x - x0, 2) + pow(z - z0, 2)) / n_abs_3, - (y - y0) * (z - z0) / n_abs_3 };
        vec2 = { -(x - x0) * (z - z0) / n_abs_3, - (y - y0) * (z - z0) / n_abs_3, (pow(x - x0, 2) + pow(y - y0, 2)) / n_abs_3 };
    } else if (shape == InitShape::plane) {
        vec0 = vec1 = vec2 = {0, 0, 0};
    } else if (shape == InitShape::paraboloid) {
        double stretchX = params.at(0);
        double stretchZ = params.at(1);
        double n_abs = sqrt(pow(stretchX * (x - x0), 2) + 1 + pow(stretchZ * (z - z0), 2));

        vec0 = {stretchX / n_abs - pow(stretchX, 3) * pow(x - x0, 2) / pow(n_abs, 3), 0, - stretchX * pow(stretchZ, 2) *(x-x0) * (z-z0)/pow(n_abs, 3)};
        vec1 = {-pow(stretchX, 2) * (x - x0) / pow(n_abs, 3), 0, -pow(stretchZ, 2) * (z - z0)/pow(n_abs, 3)};
        vec2 = {-pow(stretchX, 2) * stretchZ * (x - x0) * (z - z0)/pow(n_abs, 3), 0, stretchZ / n_abs - pow(stretchZ, 3) * pow(z - z0, 2) / pow(n_abs, 3)};
    } else if (shape == InitShape::ellipsoid) {
        double stretchX = params.at(0);
        double stretchY = params.at(1);
        double stretchZ = params.at(2);
        double n_abs = sqrt(pow(stretchX * (x - x0), 2) + pow(stretchY * (y - y0), 2) + pow(stretchZ * (z - z0), 2));

        vec0 = { stretchX / n_abs - pow(stretchX, 3) * pow(x - x0, 2) / pow(n_abs, 3),
                 stretchX * pow(stretchY, 2) * (x - x0) * (y - y0) / pow(n_abs, 3),
                 stretchX * pow(stretchZ, 2) * (x - x0) * (z - z0) / pow(n_abs, 3) };
        vec1 = { -pow(stretchX, 2) * stretchY * (x - x0) * (y - y0) / pow(n_abs, 3),
                 stretchY / n_abs - pow(stretchY, 3) * pow(y - y0, 2) / pow(n_abs, 3),
                 -stretchY * pow(stretchZ, 2) * (y - y0) * (z - z0) / pow(n_abs, 3)};
        vec2 = { -pow(stretchX, 2) * stretchZ * (x - x0) * (z - z0) / pow(n_abs, 3),
                 -pow(stretchY, 2) * stretchZ * (y - y0) * (z - z0) / pow(n_abs, 3),
                 stretchZ / n_abs - pow(stretchZ, 3) * pow(z - z0, 2) / pow(n_abs, 3)};
    }

    Matrix ret = {vec0, vec1, vec2};
    return ret;
}


/**
 * Given the contact angle, calculate the normal vector.
 *
 * This function is only applicable in 2D.
 * @param initAngle The angle in radiants.
 */
Vector LevelSet::normalVector2D(double initAngle) {
    Vector initNormal;
    if (trackedCP == "left") {
        initNormal = { -sin(initAngle), cos(initAngle), 0};
    } else {
        initNormal = {sin(initAngle), cos(initAngle), 0};
    }
    return initNormal;
}

/**
 * Given the contact angle, calculate the normal vector.
 *
 * This function is only applicable in 2D. This is the non-member variant of normalVector2D(initAngle).
 * @param initAngle The angle in radiants.
 */
Vector normalVector2D(double initAngle, std::string trackedCP) {
    Vector initNormal;
    if (trackedCP == "left") {
        initNormal = { -sin(initAngle), cos(initAngle), 0};
    } else {
        initNormal = {sin(initAngle), cos(initAngle), 0};
    }
    return initNormal;
}

/**
 * Calculate the next timestep at a given time and evolves the field.
 *
 * For the flux calculation, the upwind method is used to ensure stability. Furthermore we impose
 * the Neumann boundary condition onto the field
 *
 * @param dt The width of the timestep by which to evolve the field
 * @param timestep The index of the timestep
 */
void LevelSet::calculateNextTimestep(double dt, int timestep) {
    LevelSet tempPhi(*this);

    const Vector upNormal = {0, 1, 0};
    const Vector downNormal = {0,-1, 0};
    const Vector leftNormal = {-1, 0, 0};
    const Vector rightNormal = {1, 0, 0};
    const Vector frontNormal = {0, 0, 1};
    const Vector backNormal = {0, 0, -1};

    // loop over all cells
#pragma omp parallel shared(tempPhi)
    {
#pragma omp for collapse(3)
    for (int k = 0; k < numZ; k++) {
        for (int j = 0; j < numY; j++) {
            for (int i = 0; i < numX; i++) {

                //Calculate the flux of phi over the cell faces
                double flux = 0;
                double sp = 0;

                for (int dir = 0; dir < 6; dir++) {

                    switch(dir) {
                    case 0:
                        sp = field->at(timestep*dt, (i + 0.5)*dx, (j + 1)*dy, (k + 0.5)*dz) * upNormal;

                        if (j == numY - 1) {
                            flux += sp*tempPhi.at(i, j, k)*dx*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j + 1, k)) * dx * dz;
                        }
                        break;

                    case 1:
                        sp = field->at(timestep*dt, (i + 0.5) * dx, j * dy, (k + 0.5) * dz) * downNormal;
                        if (j == 0 ) {
                            flux += sp * tempPhi.at(i, j, k) * dx * dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j - 1, k)) * dx * dz;
                        }
                        break;

                    case 2:
                        sp = field->at(timestep*dt, i * dx, (j + 0.5) * dy, (k + 0.5) * dz) * leftNormal;
                        if (i == 0) {
                            flux += sp * tempPhi.at(i, j, k) * dy * dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i - 1, j, k)) * dy * dz;
                        }
                        break;

                    case 3:
                        sp = field->at(timestep*dt, (i + 1) * dx, (j + 0.5) * dy, (k + 0.5) * dz) * rightNormal;
                        if (i == numX - 1) {
                            flux += sp * tempPhi.at(i, j, k) * dy * dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i + 1, j, k)) * dy * dz;
                        }
                        break;

                    case 4:
                        if (numZ == 1)
                            break; // only relevant for 3D
                        sp = field->at(timestep*dt, (i + 0.5) * dx, (j + 0.5) * dy, (k + 1) * dz) * frontNormal;
                        if (k == numZ - 1) {
                            flux += sp * tempPhi.at(i, j, k) * dx * dy;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j, k + 1)) * dx * dy;
                        }
                        break;

                    case 5:
                        if (numZ == 1)
                            break; // only relevant for 3D
                        sp = field->at(timestep*dt, (i + 0.5) * dx, (j + 0.5) * dy, k * dz) * backNormal;
                        if (k == 0) {
                            flux += sp * tempPhi.at(i, j, k) * dx * dy;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j, k - 1)) * dx * dy;
                        }
                        break;
                    }
                }
                this->at(i, j, k) = this->at(i, j, k) - dt / (dx * dy * dz) * flux;
            }
        }
    }

    }
}

/**
 * Calculate the next timestep at a given time and evolves the field using
 * our new surface term approach to preserve the norm of the gradient at the
 * zero level set
 *
 * For the flux calculation, the upwind method is used to ensure stability. Furthermore we impose
 * the Neumann boundary condition onto the field
 *
 * @param dt The width of the timestep by which to evolve the field
 * @param timestep The index of the timestep
 */
void LevelSet::calculateNextTimestepSourceTerm(double dt, int timestep) {
    LevelSet tempPhi(*this);

    const Vector upNormal = {0, 1, 0};
    const Vector downNormal = {0,-1, 0};
    const Vector leftNormal = {-1, 0, 0};
    const Vector rightNormal = {1, 0, 0};
    const Vector frontNormal = {0, 0, 1};
    const Vector backNormal = {0, 0, -1};

    Vector reconstructedNormal = {0,0,0};

    // loop over all cells
#pragma omp parallel shared(tempPhi)
    {
#pragma omp for collapse(3)
        for (int k = 0; k < numZ; k++) {
            for (int j = 0; j < numY; j++) {
                for (int i = 0; i < numX; i++) {

                    //Calculate the flux of phi over the cell faces
                    double flux = 0;
                    double sp = 0;
                    double source = 0;

                    for (int dir = 0; dir < 6; dir++) {

                        switch(dir) {
                            case 0:
                                sp = field->at(timestep*dt, (i + 0.5)*dx, (j + 1)*dy, (k + 0.5)*dz) * upNormal;
                                if (j == numY - 1) {
                                    flux += sp*tempPhi.at(i, j, k)*dx*dz;
                                }
                                else{
                                    flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j + 1, k)) * dx * dz;
                                }
                                break;

                            case 1:
                                sp = field->at(timestep*dt, (i + 0.5) * dx, j * dy, (k + 0.5) * dz) * downNormal;
                                if (j == 0 ) {
                                    flux += sp * tempPhi.at(i, j, k) * dx * dz;
                                }
                                else{
                                    flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j - 1, k)) * dx * dz;
                                }
                                break;

                            case 2:
                                sp = field->at(timestep*dt, i * dx, (j + 0.5) * dy, (k + 0.5) * dz) * leftNormal;
                                if (i == 0) {
                                    flux += sp * tempPhi.at(i, j, k) * dy * dz;
                                }
                                else{
                                    flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i - 1, j, k)) * dy * dz;
                                }
                                break;

                            case 3:
                                sp = field->at(timestep*dt, (i + 1) * dx, (j + 0.5) * dy, (k + 0.5) * dz) * rightNormal;
                                if (i == numX - 1) {
                                    flux += sp * tempPhi.at(i, j, k) * dy * dz;
                                }
                                else{
                                    flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i + 1, j, k)) * dy * dz;
                                }
                                break;

                            case 4:
                                if (numZ == 1)
                                    break; // only relevant for 3D
                                sp = field->at(timestep*dt, (i + 0.5) * dx, (j + 0.5) * dy, (k + 1) * dz) * frontNormal;
                                if (k == numZ - 1) {
                                    flux += sp * tempPhi.at(i, j, k) * dx * dy;
                                }
                                else{
                                    flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j, k + 1)) * dx * dy;
                                }
                                break;

                            case 5:
                                if (numZ == 1)
                                    break; // only relevant for 3D
                                sp = field->at(timestep*dt, (i + 0.5) * dx, (j + 0.5) * dy, k * dz) * backNormal;
                                if (k == 0) {
                                    flux += sp * tempPhi.at(i, j, k) * dx * dy;
                                }
                                else{
                                    flux += (fmax(sp,0.0)*tempPhi.at(i, j, k) + fmin(sp, 0.0) * tempPhi.at(i, j, k - 1)) * dx * dy;
                                }
                                break;
                        }
                    }

                    // Compute the local normal (2D)
                    if(i==0){
                      reconstructedNormal[0] = (tempPhi.at(i+1,j,k)-tempPhi.at(i,j,k))/(dx);
                    }
                    else if(i==numX-1){
                      reconstructedNormal[0] = (tempPhi.at(i,j,k)-tempPhi.at(i-1,j,k))/(dx);
                    }
                    else {
                      reconstructedNormal[0] = (tempPhi.at(i+1,j,k)-tempPhi.at(i-1,j,k))/(2*dx);
                    }

                    if(j==0){
                      reconstructedNormal[1] = (tempPhi.at(i,j+1,k)-tempPhi.at(i,j,k))/(dy);
                    }
                    else if(j==numY-1){
                      reconstructedNormal[1] = (tempPhi.at(i,j,k)-tempPhi.at(i,j-1,k))/(dy);
                    }
                    else {
                      reconstructedNormal[1] = (tempPhi.at(i,j+1,k)-tempPhi.at(i,j-1,k))/(2*dy);
                    }

                    reconstructedNormal[2] = 0;

                    reconstructedNormal = reconstructedNormal/abs(reconstructedNormal);


                    // Compute source term
                    source = (field->gradAt(timestep*dt, (i + 0.5)*dx, (j + 0.5)*dy, (k + 0.5)*dz)*reconstructedNormal)*reconstructedNormal*(-1.0);

                    // explicit update
                    this->at(i, j, k) = this->at(i, j, k)*(1.0-source*dt) - dt / (dx * dy * dz) * flux;

                    // implicit update
                    //this->at(i, j, k) = (this->at(i, j, k) - dt / (dx * dy * dz) * flux)/(1+source*dt);

                }
            }
        }

    }
}

double LevelSet::getGradPhiNormAtContactPoint(int timestep) {
    double minimum;
    std::array<int, 3> contactPoint = getContactPointIndices(timestep);
    Vector gradPhi = getNormalVector(contactPoint, false, false);
    double gradPhiNorm = abs(gradPhi);
    return gradPhiNorm;
}


std::vector<Vector> LevelSet::getPositionReference() const {
    return positionReference;
}

std::vector<double> LevelSet::getAngleReference() const {
    return angleReference;
}

const std::vector<double>& LevelSet::getCurvatureReference() const {
    return curvatureReference;
}

std::vector<Vector> LevelSet::getTangentAReference() const {
    return tangentAReference;
}

std::vector<Vector> LevelSet::getTangentBReference() const {
    return tangentBReference;
}

std::vector<Vector> LevelSet::getTangentCReference() const {
    return tangentCReference;
}

std::vector<double> LevelSet::getSectionalCurvatureAReference() const {
    return sectionalCurvatureAReference;
}

std::vector<double> LevelSet::getSectionalCurvatureBReference() const {
    return sectionalCurvatureBReference;
}

std::vector<double> LevelSet::getSectionalCurvatureCReference() const {
    return sectionalCurvatureCReference;
}
