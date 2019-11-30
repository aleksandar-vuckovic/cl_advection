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
        std::string trackedCP, double dt, int timesteps, array<double, 3> expCP, array<double, 3> expNormalVec, double initCurvature)
        : Field<double>(numX, numY, numZ), positionReference(timesteps), angleReference(timesteps), curvatureReference(timesteps){
		this->dx = dx;
		this->dy = dy;
		this->dz = dz;
		this->field = field;
		this->trackedCP = trackedCP;


	    array<double, 3> initCP = getInitCP(expCP, 0.001);
        double expAngle = acos(expNormalVec[1]);
		// Calculate reference data
        contactPointExplicitEuler(dt, timesteps, initCP);
        if (numZ == 1 && (field->getName() == "navierField" || field->getName() == "timeDependentNavierField")) {
            referenceAngleLinearField(dt, timesteps, expAngle);
        } else {
            referenceAngleExplicitEuler(dt, timesteps, expNormalVec);
        }

        if (field->getName() == "quadraticField") {
            referenceCurvatureQuadraticField(dt, timesteps, initCurvature);
        } else {
            referenceCurvatureExplicitEuler(dt, timesteps, initCurvature,  expAngle);
        }
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
array<double, 3> LevelSet::getInitCP(array<double, 3> expcp, double epsilon) {
    array<double, 3> candidate = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
        for (int y = 0; y < this->numY; y++)
            for (int z = 0; z < this->numZ; z++) {
                array<double, 3> other = {x*dx, y*dy, z*dz};
                if (abs(candidate - expcp) > abs(other - expcp) && std::abs(this->at(x, y, z)) < epsilon && y == 0) {
                    candidate = other;
		}
	}
    return candidate;
}

/**
 * Calculate the new position of the contact point at a given time.
 * Using the explicit euler method, calculate the position of the contact point initCP at all times and writes them to
 * the reference array.
 *
 * @param dt The length of a single timestep
 * @param timestep The index of the timestep
 * @param initCP The position of the contact point at timestep 0
 */
void LevelSet::contactPointExplicitEuler(double dt, int last_timestep, array<double, 3> initCP) {
    array<double, 3> &temp = initCP;
    for (int i = 0; i < last_timestep; i++) {
    	temp = temp + dt*field->at(i*dt, temp[0], temp[1], temp[2]);
    	positionReference[i] = temp;
    }
}

/**
 * Calculate the contact point position for the 15navier field using the analytic solution.
 *
 * @param t The time at which to calculate the contact point position
 * @param c1 A parameter of the navier field
 * @param x0 The initial position of the contact point
 * @param v0 A parameter of the navier field
 * @return The position of the contact point in arbitrary units
 */
array<double, 3> LevelSet::contactPointLinearField(double t, double c1, double x0, double v0) {
	if (field->getName() == "navierField") {
		return {x0 * exp(c1 * t) + v0/c1 * (exp(c1 * t) - 1), 0, 0};
	} else if (field->getName() == "timeDependentNavierField") {
		double tau = field->getTau();
		return {x0 *exp(c1/M_PI * tau*sin(M_PI*t/tau)) + v0/c1 * (exp(c1/M_PI * tau*sin(M_PI*t/tau)) -1), 0, 0};
	} else {
		throw std::invalid_argument("Please choose either navierField or timeDependentNavierField if analyzing linear fields.");
	}
}

/**
 * Return the indices of the contact point in 2D or the coordinates most closely matching the given parameter.
 * This function works differently for 2D and 3D. If the simulation is 2D, it ignores the input parameter and returns
 * the contact point by checking where the sign of the LevelSet field changes.
 *
 * In 3D, it returns the coordinates of the point according to the reference ODE. In other words, it is assumed that in 3D,
 * the solution exactly matches the actual contact point. It uses a convex combination of 
 * neighboring LevelSet values depending on whether the right or the left contact point is being tracked.
 *
 * @param point A point in the space of the Level set field.
 * @param indexOnly Whether to return only the indices of the point or its actual coordinates. Default is false.
 * @return The contact point
 */
array<double, 3> LevelSet::getContactPoint(int timestep, bool indexOnly /* = false */) {
    if (numZ == 1) {
        if (trackedCP == "left") {
            double initSign = this->at(0, 0, 0)/std::abs(this->at(0, 0, 0));
            for (int i = 0; i < numX; i++)
                if (this->at(i, 0, 0)*initSign < 0) {
                    
                    if (indexOnly)
                        return {double(i) + 0.5, 0, 0}; //Add 0.5 to compensate floating-point errors
                    double alpha = this->at(i,0,0)-this->at(i -1,0,0);
                    if(std::abs(alpha) < 1E-12) {
                        throw std::runtime_error("LevelSet::getContactPoint: \nDifference of LevelSet values at contact point too small for convex combination");
                    }

                    alpha = this->at(i,0,0)/(alpha);
                    return {(1 - alpha)*i*dx + alpha*(i - 1)*dx, 0, 0};
                }
        } else {
            double initSign = this->at(numX - 1, 0, 0)/std::abs(this->at(numX - 1, 0, 0));
            for (int i = numX - 1; i >= 0; i--) {
                if (this->at(i, 0, 0)*initSign < 0) {

                    if (indexOnly)
                        return {double(i) + 0.5, 0, 0};
                    double alpha = this->at(i + 1, 0, 0) - this->at(i, 0, 0);
                    if(std::abs(alpha) < 1E-12) {
                        throw std::runtime_error("LevelSet::getContactPoint: \nDifference of LevelSet values at contact point too small for convex combination");
                    }

                alpha = this->at(i + 1, 0, 0) / alpha;
                return {alpha*i*dx + (1 - alpha)*(i+1)*dx, 0, 0};
                }
            }
        }
    } else {
        if (indexOnly) {
            array<double, 3> cell = {0, 0, 0};
            const array<double, 3> ref = positionReference[timestep];
            for (int i = 0; i < numX; ++i)
                for (int k = 0; k < numZ; ++k) {
                    array<double, 3> current = {i*dx, 0, k*dz};
                    array<double, 3> temp = {int(cell[0]) * dx, 0, int(cell[2]) * dz};
                    if ( abs(current - ref) < abs(temp - ref))
                        cell = {double(i) + 0.5, 0, double(k) + 0.5};
                }
            return cell;
        } else {
            return positionReference[timestep];
        }
    }
    throw std::runtime_error("LevelSet::getContactPoint: \nReached end of function without returning a value");
}

/**
 * Get the current indices of the contact point.
 * 
 * This is a wrapper function for LevelSet::getContactPoint.
 * @param timestep The timestep
 * @return The indices of the contact point.
 */
array<int, 3> LevelSet::getContactPointIndices(int timestep) {
    array<double, 3> vec = getContactPoint(timestep, true);
    return { int(vec[0]), int(vec[1]), int(vec[2]) };
}

/**
 * Calculate the normal vector.
 * 
 * This function uses finite differences to calculate the normal vector of the Level set field at cell
 * @param cell The indices of the contact point.
 * @return The normal vector.
 */
array<double, 3> LevelSet::getNormalVector(array<int, 3> cell) const {
    double normalX = 0, normalY = 0, normalZ = 0;

    if (numZ == 1) {

	if (trackedCP == "left") {
		 // find root of phi, alpha: coefficient for convex combination
		 double alpha = this->at(cell[0],0,0)-this->at(cell[0]-1,0,0);

		 if(std::abs(alpha) < 1E-12){
		   throw std::runtime_error("LevelSet::getContactAngle: \nDifference of LevelSet values at contact point too small for convex combination");
		 }

		 alpha = this->at(cell[0],0,0)/(alpha);

		//Calculate angle at this cell with finite differences
		normalX = alpha*(this->at(cell[0], cell[1], cell[2]) - this->at(cell[0]-2, cell[1], cell[2]))/(2*dx)
		+ (1-alpha)*(this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2]))/(2*dx);
		normalY = alpha*(-this->at(cell[0]-1, cell[1] + 2, cell[2])
						  + 4.0*this->at(cell[0]-1, cell[1] + 1, cell[2])
						  - 3.0*this->at(cell[0]-1, cell[1], cell[2]))/(2*dy)
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
                normalX = (this->at(cell[0] + 1, cell[1], cell[2]) - this->at(cell[0] - 1, cell[1], cell[2]))/(2*dx);
        } else if (cell[0] == 0) {
                normalX = (-1*this->at(cell[0] + 2, cell[1], cell[2]) +4*this->at(cell[0] + 1, cell[1], cell[2]) -3*this->at(cell[0], cell[1], cell[2]))/(2*dx);
        } else {
                normalX = (this->at(cell[0] - 2, cell[1], cell[2]) -4*this->at(cell[0] - 1, cell[1], cell[2]) +3*this->at(cell[0], cell[1], cell[2]))/(2*dx);
        }

        normalY = (-1*this->at(cell[0], cell[1] + 2, cell[2]) +4*this->at(cell[0], cell[1] + 1, cell[2]) -3*this->at(cell[0], cell[1], cell[2]))/(2*dy);

        if (cell[2] > 0 && cell[2] < numZ - 1) {
                normalZ = (this->at(cell[0], cell[1], cell[2] + 1) - this->at(cell[0], cell[1], cell[2] - 1))/(2*dz);
        } else if (cell[2] == 0) {
                normalZ = (-1*this->at(cell[0], cell[1], cell[2] + 2) +4*this->at(cell[0], cell[1] + 1, cell[2]) -3*this->at(cell[0], cell[1], cell[2]))/(2*dz);
        } else {
                normalZ = (this->at(cell[0], cell[1], cell[2] - 2) -4*this->at(cell[0], cell[1], cell[2] - 1) +3*this->at(cell[0], cell[1], cell[2]))/(2*dz);
        }
    }

    array<double ,3> normal = {normalX, normalY, normalZ};
    normal = normal/abs(normal);
    
    return normal;
}

/**
 * Calculate the contact angle.
 * This function uses finite differences to calculate the normal vector of the Level set field at cell and
 * thus calculate the contact angle. This is a wrapper function for LevelSet::getNormalVector. 
 *
 * @param cell The indices of the contact point.
 * @return The contact angle in degrees
 */
double LevelSet::getContactAngle(array<int, 3> cell) {
    array<double, 3> normal = getNormalVector(cell);
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
void LevelSet::referenceAngleExplicitEuler(double dt, int last_timestep, array<double, 3> n_sigma_init) {
	array<double, 3> &n_sigma = n_sigma_init;
	array<double, 3> deriv = {0, 0, 0};
	for (int i = 0; i < last_timestep; i++) {
	    angleReference[i] = acos(n_sigma[1])/M_PI*180;
	    array<double, 3> CP = positionReference[i];
		deriv = -1*transpose(field->gradAt(i*dt, CP[0], CP[1], CP[2]))*n_sigma + ((field->gradAt(i*dt, CP[0], CP[1], CP[2])*n_sigma)*n_sigma)*n_sigma;
		n_sigma = deriv*dt + n_sigma;
		n_sigma = n_sigma/abs(n_sigma);
	}
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
 *
 * @param dt The length of a timestep
 * @param timestep The index of the timestep.
 * @param initCurvature The initial curvature at time t == 0
 * @param CP The coordinates of the contact point
 * @param cell The indices of the contact point
 */
void LevelSet::referenceCurvatureExplicitEuler(double dt, int last_timestep, double initCurvature, double initAngle) {
    double curvature = initCurvature;
    array<double, 3> initNormal;
    if (trackedCP == "left") {
        initNormal = { -sin(initAngle), cos(initAngle)};
    } else {
        initNormal = {sin(initAngle), cos(initAngle)};
    }

    for (int i = 0; i < last_timestep; i++) {
        curvatureReference[i] = curvature;
        array<double, 3> CP = positionReference[i];
        double contactAngle = angleReference[i]/180*M_PI;
        array<double, 3> normal, tau;
        if (trackedCP == "left") {
            normal = { -sin(contactAngle), cos(contactAngle), 0};
            tau = {cos(contactAngle), sin(contactAngle), 0};
        } else {
            normal = {sin(contactAngle), cos(contactAngle), 0};
            tau = { -cos(contactAngle), sin(contactAngle), 0};
        }

        // temp is the second derivative of v in the tau direction (= y direction)
        array<double, 3> temp;
        if (field->getName() == "shearField") {
            double v0 = field->getV0();
            array<double, 3> row1 = {v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1])*tau[0] + v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1])*tau[1],
                                    v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1])*tau[0] + v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1])*tau[1],
                                    0};
            array<double, 3> row2 = {-v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1])*tau[0] - v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1])*tau[1],
                                     -v0*M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1])*tau[0]   -v0*M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1])*tau[1]};
            array<double, 3> row3 = {0, 0, 0};
            array<array<double, 3>, 3> M = {row1, row2, row3};
            temp = M*tau;
        } else if (field->getName() == "navierField") {
            temp = {0, 0, 0};
        }

        curvature = curvature + dt*(temp*normal - 3*curvature*(field->gradAt(i*dt, CP[0], CP[1], CP[2])*tau)*tau);
    }
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
double LevelSet::getCurvatureDivergence(array<int, 3> cell) const {
    /* Define what number of cells in each direction(excluding the main cell) are considered local.
       The resulting normal vector field is defined on a cuboid with sidelength 2*local + 1 */
    int local = 2;
    int sidelengthZ;
    if (numZ > 1)
        sidelengthZ = 2*local + 1;
    else
        sidelengthZ = 1;

    // Declare a field of normal vectors
    Field<array<double, 3> > localField(2*local + 1, local + 1, sidelengthZ);

    for (int x = -local; x <= local; x++)
        for (int y = 0; y <= local; y++)
            for (int z = sidelengthZ/2; z <= sidelengthZ/2; z++) {
                array<int, 3> temp = {x, y, z};
                temp = temp + cell;

                array<double, 3> normal = getNormalVector(cell);
                if (this->numZ > 1)
                    localField.at(x+local, y, z+local) = normal;
                else
                    localField.at(x+local, y, 0) = normal;
            }

    // Calculate divergence of localField at cell, which by definition
    // is in the "middle" of localField
    double dnx_dx = (localField.at(local + 1, 0, sidelengthZ/2)[0] - localField.at(local - 1, 0, sidelengthZ/2)[0]) / (2*dx);

    double dnx_dy = (-localField.at(local, 2, sidelengthZ/2)[0]
                + 4.0*localField.at(local, 1, sidelengthZ/2)[0]
                - 3.0*localField.at(local, 0, sidelengthZ/2)[0] ) / (2*dy);

    double dny_dx = (localField.at(local + 1, 0, sidelengthZ/2)[1] - localField.at(local - 1, 0, sidelengthZ/2)[1]) / (2*dx);

    double dny_dy = (-localField.at(local, 2, sidelengthZ/2)[1]
                + 4.0*localField.at(local, 1, sidelengthZ/2)[1]
                - 3.0*localField.at(local, 0, sidelengthZ/2)[1] ) / (2*dy);

    double dnx_dz, dny_dz, dnz_dx, dnz_dy, dnz_dz;
    dnx_dz = dny_dz = dnz_dx = dnz_dy = dnz_dz = 0;

    if (numZ > 1) {
        dnx_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[0] - localField.at(local, 0, sidelengthZ/2 - 1)[0]) / (2*dz);

        dny_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[1] - localField.at(local, 0, sidelengthZ/2 - 1)[1]) / (2*dz);

        dnz_dx = (localField.at(local + 1, 0, sidelengthZ/2)[2] - localField.at(local - 1, 0, sidelengthZ/2)[2]) / (2*dx);

        dnz_dy = (-localField.at(local, 2, sidelengthZ/2)[2]
                    + 4.0*localField.at(local, 1, sidelengthZ/2)[2]
                    - 3.0*localField.at(local, 0, sidelengthZ/2)[2] ) / (2*dy);

        dnz_dz = (localField.at(local, 0, sidelengthZ/2 + 1)[2] - localField.at(local, 0, sidelengthZ/2 - 1)[2]) / (2*dz);
    }
    
    array<double, 3> normal = localField.at(local, 0, sidelengthZ/2);

    // This still only works for 2D
    array<double, 3> tau = {normal[1], -normal[0], 0};

    array<double, 3> row1 = {dnx_dx, dnx_dy, dnx_dz};
    array<double, 3> row2 = {dny_dx, dny_dy, dny_dz};
    array<double, 3> row3 = {dnz_dx, dnz_dy, dnz_dz};

    array< array<double, 3>, 3> gradNormal = {row1, row2, row3};

    double kappa = -1*(gradNormal*tau)*tau;

    return kappa;
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
 * @param xmfFile A pointer to the XMF file.
 */
void LevelSet::writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile) {
    int Npoints = numX*numY*numZ;
    double *pointCoordinates = new double[Npoints*3];
    double *pointPhiValues = new double [Npoints];

    int index = 0;
    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++)
            for (int i = 0; i < numX; i++) {
                if (timestep == 0) {
                    pointCoordinates[index] = i*dx;
                    pointCoordinates[index +1] = j*dy;
                    pointCoordinates[index +2] = k*dz;
                    index += 3;
                }
                pointPhiValues[i + j*numX + k*numX*numY] = this->at(i, j, k);
	    }

    // If it is the first iteration, create coordinate file
    if (timestep == 0) {
        *xmfFile << "<?xml version=\"1.0\" ?>\n"
                 << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
                 << "<!ENTITY Npoints \"" + std::to_string(Npoints) + "\">\n"
                 << "]>\n"
                 << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
                 << "<Domain>\n"
                 << "<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
                 << "<Time TimeType=\"List\">\n"
                 << "<DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\""+ std::to_string(total_writesteps)+"\">\n";
        for (int i = 0; i < total_writesteps; i++) {
            *xmfFile << std::to_string(i*((double)total_timesteps/total_writesteps)*dt) << " ";
        }
        *xmfFile <<"</DataItem>\n" << "</Time>\n";

        //Write field coordinates into binary file
        FILE *fieldFile;
        fieldFile = fopen("data/field.bin", "wb");
        fwrite(pointCoordinates, sizeof(double), Npoints*3, fieldFile);
        fclose(fieldFile);
    }
    FILE *PhiFile;
    std::string filename = "data/Phi_t="+ std::to_string(timestep*dt)+".bin";
    PhiFile = fopen(filename.data(), "wb");
	fwrite(pointPhiValues, sizeof(double), Npoints, PhiFile);
	fclose(PhiFile);



	*xmfFile << "<Grid>\n"
			 << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\""+std::to_string(Npoints) +"\"/>\n"
			 << "<Geometry GeometryType=\"XYZ\"> \n"
			 << "<DataItem Name=\"points\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
			 << "field.bin\n"
			 << "</DataItem>\n</Geometry>\n"
			 << "<Attribute Name =\"lvlset\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
			 << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\"  Endian=\"Little\" Dimensions=\"&Npoints;\">\n"
			 << "Phi_t=" + std::to_string(timestep*dt) +".bin\n"
			 << "</DataItem></Attribute>\n"
			 << "<Attribute Name =\"VelocityField\" AttributeType=\"Vector\" Center=\"Cell\">\n"
			 << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n";

	if (field->getName() == "timeDependentNavierField")
		*xmfFile << "Vel_t=" + std::to_string(timestep*dt) +".bin\n";
	else
		*xmfFile << "Vel_t=" + std::to_string(0*dt) +".bin\n";

	*xmfFile << "</DataItem></Attribute>\n"
			 /*<< "<Attribute Name =\"Tangential Vector\" AttributeType=\"Vector\" Center=\"Cell\">\n"
             << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
             << "Tau_t=" + std::to_string(timestep*dt) +".bin\n"
             << "</DataItem></Attribute>\n"  */
			 << "<Attribute Name =\"Streamlines\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
			 << "<DataItem Format=\"Binary\" NumberType=\"Int\" Precision=\"4\" Endian=\"Little\" Dimensions=\"&Npoints;\">\n"
			 << "stream.bin\n"
			 << "</DataItem></Attribute>\n"
			 << "</Grid>\n";


    delete[] pointCoordinates;
    delete[] pointPhiValues;

    //Write velocity field
    if (field->getName() == "timeDependentNavierField")
    	field->writeToFile(timestep*dt);
    else if (timestep == 0)
    	field->writeToFile(0*dt);

    //Write tangential vector to file
   //writeTangentialVectorToFile(dt, timestep);
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
    double *fieldValues = new double[numX*numY*numZ*3];
    int index = 0;

    array<int, 3> cell = getContactPointIndices(timestep);
    
    
    for (int k = 0; k < numZ; k++) {
        for (int j = 0; j < numY; j++) {
            for (int i = 0; i < numX; i++) {
                array<double, 3> temp;
                if (i == cell[0] && j == cell[1] && k == cell[2]) {
                    temp = getNormalVector(cell);
                } else {
                    temp = {0, 0, 0};
                }
                fieldValues[index] = temp[0];
                fieldValues[index + 1] = temp[1];
                fieldValues[index + 2] = temp[2];
                index += 3;
            }
        }
    }
    
    
    std::string filename = "data/Tau_t=" + std::to_string(dt*timestep) + ".bin";
    FILE *tauFile;
    tauFile = fopen(filename.data(), "wb");
    fwrite(fieldValues, sizeof(double), numX*numY*numZ*3, tauFile);
    fclose(tauFile);

    delete[] fieldValues;
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
void LevelSet::initDroplet(array<double, 3> center, double radius) {
    for (int x = 0; x < this->numX; x++)
        for (int y = 0; y < this->numY; y++)
            for (int z = 0; z < this->numZ; z++)
                this->at(x, y, z) = pow(x*dx - center[0], 2) + pow(y*dy - center[1], 2) + pow(z*dz - center[2], 2) - pow(radius, 2);
}

/** Initialize a plane.
 * 
 *  @param refPoint The reference point of the plane in degrees.
 *  @param angleA The angle between the initialized plane and the x-z-plane (= polar angle of normal vector).
 *  @param angleB The angle of rotation around the y-axis (= azimuthal angle).
 */
void LevelSet::initPlane(array<double, 3> refPoint, double polarAngle, double azimuthalAngle) {
    polarAngle = polarAngle/180*M_PI;
    azimuthalAngle = azimuthalAngle/180*M_PI;
    array<double, 3> normal = {sin(polarAngle) * cos(azimuthalAngle), cos(polarAngle), sin(polarAngle)* sin(azimuthalAngle)};
    
    for (int x = 0; x < numX; x++)
        for (int y = 0; y < numY; y++)
            for (int z = 0; z < numZ; z++) 
                this->at(x, y, z) = normal * (array<double, 3>({x*dx, y*dy, z*dz}) - refPoint);
            
}

/**
 * Given the contact angle, calculate the normal vector.
 *
 * This function is only applicable in 2D.
 * @param initAngle The angle in radiants.
 */
array<double, 3> LevelSet::normalVector2D(double initAngle) {
    array<double, 3> initNormal;
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
array<double, 3> normalVector2D(double initAngle, std::string trackedCP) {
    array<double, 3> initNormal;
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

    const array<double, 3> upNormal = {0, 1, 0};
    const array<double, 3> downNormal = {0,-1, 0};
    const array<double, 3> leftNormal = {-1, 0, 0};
    const array<double, 3> rightNormal = {1, 0, 0};
    const array<double, 3> frontNormal = {0, 0, 1};
    const array<double, 3> backNormal = {0, 0, -1};

    // loop over all cells
#pragma omp parallel shared(tempPhi)
    {
#pragma omp for collapse(3)
    for (int x = 0; x < numX; x++) {
        for (int y = 0; y < numY; y++) {
            for (int z = 0; z < numZ; z++) {

                //Calculate the flux of phi over the cell faces
                double flux = 0;
                double sp = 0;

                for (int dir = 0; dir < 6; dir++) {

                    switch(dir) {
                    case 0:
                        sp = field->at(timestep*dt, (x+1/2)*dx, (y+1)*dy, (z+1/2)*dz) * upNormal;
                        if (y == numY - 1) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y+1, z))*dx*dz;
                        }
                        break;

                    case 1:
                        sp = field->at(timestep*dt, (x+1/2)*dx, y*dy, (z+1/2)*dz) * downNormal;
                        if (y == 0 ) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y-1, z))*dx*dz;
                        }
                        break;

                    case 2:
                        sp = field->at(timestep*dt, x*dx, (y+1/2)*dy, (z+1/2)*dz) * leftNormal;
                        if (x == 0) {
                            flux += sp*tempPhi.at(x, y, z)*dy*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x-1, y, z))*dy*dz;
                        }
                        break;

                    case 3:
                        sp = field->at(timestep*dt, (x+1)*dx, (y+1/2)*dy, (z+1/2)*dz) * rightNormal;
                        if (x == numX - 1) {
                            flux += sp*tempPhi.at(x, y, z)*dy*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x+1, y, z))*dy*dz;
                        }
                        break;

                    case 4:
                        if (numZ == 1)
                            break; // only relevant for 3D
                        sp = field->at(timestep*dt, (x+1/2)*dx, (y+1/2)*dy, (z+1)*dz) * frontNormal;
                        if (z == numZ - 1) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dy;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y, z+1))*dx*dy;
                        }
                        break;

                    case 5:
                        if (numZ == 1)
                            break; // only relevant for 3D
                        sp = field->at(timestep*dt, (x+1/2)*dx, (y+1/2)*dy, z*dz) * backNormal;
                        if (z == 0) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dy;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y, z-1))*dx*dy;
                        }
                        break;
                    }
                }
                this->at(x, y, z) = this->at(x, y, z) - dt/(dx*dy*dz)*flux;
            }
        }
    }

    }
}

std::vector<array<double, 3>> LevelSet::getPositionReference() {
    return positionReference;
}

std::vector<double> LevelSet::getAngleReference() {
    return angleReference;
}

std::vector<double> LevelSet::getCurvatureReference() {
    return curvatureReference;
}
