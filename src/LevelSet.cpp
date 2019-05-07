#include "LevelSet.hpp"

array<double, 3> LevelSet::getInitCP(double dt, array<double, 3> expcp, double epsilon) {
    array<double, 3> candidate = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++) {
		array<double, 3> other = {x*dx, y*dx, z*dx};
		if (abs(candidate - expcp) > abs(other - expcp) && std::abs(this->at(x, y, z)) < epsilon && y == 0) {
		    candidate = other;
		}
	    }

    return candidate;
}

array<double, 3> LevelSet::getContactPoint(double dt, int timestep, int timesteps, array<double, 3> initCP) {
    array<double, 3> &temp = initCP;
    //Calculate current position of contact point
    for (int i = 0; i < timestep; i++)
	temp = temp + dt*field->at(temp[0], temp[1], temp[2]);

    return temp;
}

array<int, 3> LevelSet::getContactPointCoordinates(array<double, 3> point) {
    // If simulation is 2D, ignore input parameter and return
    // coordinates of /left/ contact point by checking where the sign of the
    // LevelSet field changes
    if (this->numZ == 1){
        double initSign = this->at(0, 0, 0)/abs(this->at(0, 0, 0));
        for (int x = 0; x < this->numX; x++)
            if (this->at(x, 0, 0)*initSign < 0)
                return {x, 0, 0};
    } else {
        //Find cell corresponding to this point
        array<int, 3> cell = {0, 0, 0};
        for (int x = 0; x < this->numX; x++)
            for (int y = 0; y < this->numY; y++)
                for (int z = 0; z < this->numZ; z++) {
                    array<int, 3> other = {x, y, z};
                    if (abs(point - cell*dx) > abs(point - other*dx))
                        cell = other;
                }
        return cell;
    }
    return array<int, 3> {0,0,0};
}
     
double LevelSet::getContactAngle(double dt, double timestep, array<int, 3> cell) {
    //Calculate angle at this cell with finite differences
    double normalX = (this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2]))/(2*dx);
    double normalY = (-this->at(cell[0], cell[1] + 2, cell[2])
                      + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
                      - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dx); // second order difference quotient
    double normalZ;
    if (this->numZ > 1)
	normalZ = (this->at(cell[0], cell[1], cell[2]+1) - this->at(cell[0]-1, cell[1], cell[2]-1))/(2*dx);
    else
	normalZ = 0;
    array<double ,3> normal = {normalX, normalY, normalZ};
    normal = normal/abs(normal);
    return acos(normal[1]);
}

double LevelSet::getReferenceAngleLinearField(double t, double c1, double c2, double theta0) {
    return M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*t) + c2 * (exp(2*c1*t) - 1)/(2*c1));
}

double LevelSet::getReferenceCurvature(double dt, double timestep, double initCurvature, array<double, 3> CP, array<int, 3> cell) {
    double curvature = initCurvature;
    for (int i = 0; i < timestep; i++){
        
        //Calculate angle at this cell with finite differences
        double normalX = (this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2])) / (2*dx);
        double normalY = (-this->at(cell[0], cell[1] + 2, cell[2]) + 4.0*this->at(cell[0], cell[1] + 1, cell[2]) - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dx);
        double normalZ;
        if (this->numZ > 1)
            normalZ = (this->at(cell[0], cell[1], cell[2]+1) - this->at(cell[0]-1, cell[1], cell[2]-1))/(2*dx);
        else
            normalZ = 0;
        array<double, 3> normal = {normalX, normalY, normalZ};
        normal = normal/abs(normal);

        // This is the second derivative of v in the tau direction (= y direction)
        array<double, 3> temp = {M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1]), -M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1]), 0};
        
        array<double,3> tau = {normal[1], -normal[0], 0};
        
        curvature = curvature + dt*(temp*normal - 2*curvature*(field->gradAt(CP[0], CP[1], CP[2])*tau)*tau);
    }

    return curvature;
}

double LevelSet::getCurvature(double dt, int timestep, array<int, 3> cell) const {
    /* Define what number of cells in each direction(excluding the main cell) are considered local.
       The resulting normal vector field is defined on a cuboid with  sidelength 2*local + 1 */
    int local = 2;
    int sidelength = 2*local + 1;
    int sidelengthZ;
    if (this->numZ > 1)
        sidelengthZ = sidelength;
    else
        sidelengthZ = 1;
    // Declare a field of normal vectors
    Field<array<double, 3> > localField(sidelength, sidelength, sidelengthZ);

    for (int x = -sidelength/2; x <= sidelength/2; x++)
        for (int y = 0; y <= sidelength/2; y++)
            for (int z = -sidelength/2; z <= sidelength/2; z++) {
                if (this->numZ == 1 && z != 0) {
                    continue;
                }
                array<int, 3> temp = {x, y, z};
                temp = temp + cell;

                double normalX = (this->at(temp[0] + 1, temp[1], temp[2]) - this->at(temp[0]-1, temp[1], temp[2])) / (2*dx);
                // second order difference quotient
                double normalY = (-this->at(temp[0], temp[1] + 2, temp[2]) + 4.0*this->at(temp[0], temp[1] + 1, temp[2]) - 3.0*this->at(temp[0], temp[1], temp[2])) / (2*dx);
                double normalZ;
                if (this->numZ > 1)
                    normalZ = (this->at(temp[0], temp[1], temp[2] + 1) - this->at(temp[0], temp[1], temp[2] - 1)) / (2*dx);
                else
                    normalZ = 0;
                array<double ,3> normal = {normalX, normalY, normalZ};
                normal = normal/abs(normal);
                if (this->numZ > 1)
                    localField.at(x+local, y+local, z+local) = normal;
                else
                    localField.at(x+local, y+local, 0) = normal;
            }

    // Calculate divergence of localfield at cell, which by definition
    // is in the "middle" of localField
    double dnx_dx = (localField.at(sidelength/2 + 1, sidelength/2, sidelengthZ/2)[0] - localField.at(sidelength/2 - 1, sidelength/2, sidelengthZ/2)[0]) / (2*dx);

    double dnx_dy = (-localField.at(sidelength/2, sidelength/2 + 2, sidelengthZ/2)[0]
                          + 4.0*localField.at(sidelength/2, sidelength/2 + 1, sidelengthZ/2)[0]
                          - 3.0*localField.at(sidelength/2, sidelength/2, sidelengthZ/2)[0] ) / (2*dx);

    
    double dny_dx = (localField.at(sidelength/2 + 1, sidelength/2, sidelengthZ/2)[1] - localField.at(sidelength/2 - 1, sidelength/2, sidelengthZ/2)[1]) / (2*dx);
    
    double dny_dy = (-localField.at(sidelength/2, sidelength/2 + 2, sidelengthZ/2)[1]
                          + 4.0*localField.at(sidelength/2, sidelength/2 + 1, sidelengthZ/2)[1]
                          - 3.0*localField.at(sidelength/2, sidelength/2, sidelengthZ/2)[1] ) / (2*dx);
    //2D only
    array<double, 3> normal = localField.at(sidelength/2, sidelength/2, sidelengthZ/2);
    array<double, 3> tau = {normal[1], -normal[0], 0};
    array<double, 3> row1 = {dnx_dx, dnx_dy, 0};
    array<double, 3> row2 = {dny_dx, dny_dy, 0};
    array<double, 3> row3 = {0,      0,      0};
    //For this line, eclipse is complaining even though it compiles just fine
    array< array<double, 3>, 3> gradNormal = {row1, row2, row3}; // @suppress("Invalid arguments")

    double kappa = -1*(gradNormal*tau)*tau;

    return kappa;
}



void LevelSet::writeToFile(double epsilon, double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile) {
    int Npoints = numX*numY*numZ;
    double *pointCoordinates = new double[Npoints*3];
    double *pointPhiValues = new double [Npoints];

    int index = 0;
    for (int x = 0; x < numX; x++)
        for (int y = 0; y < numY; y++)
            for (int z = 0; z < numZ; z++) {
                if (timestep == 0) {
                    pointCoordinates[index] = z*dx;
                    pointCoordinates[index +1] = y*dx;
                    pointCoordinates[index +2] = x*dx;
                    index += 3;
                }
                pointPhiValues[x + y*numX + z*numX*numY] = this->at(x, y, z);
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
            *xmfFile << std::to_string(i*(total_timesteps/total_writesteps)*dt) << " ";
        }
        *xmfFile <<"</DataItem>\n" << "</Time>\n";

        //Write field coordinates into binary file
        FILE *fieldFile;
        fieldFile = fopen("data/field.bin", "wb");
        fwrite(pointCoordinates, sizeof(double), Npoints*3, fieldFile);
    }
    FILE *PhiFile;
    std::string filename = "data/Phi_t="+ std::to_string(timestep*dt)+".bin";
    PhiFile = fopen(filename.data(), "wb");
    fwrite(pointPhiValues, sizeof(double), Npoints, PhiFile);

    *xmfFile << "<Grid>\n"
             << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\""+std::to_string(Npoints) +"\"/>\n"
             << "<Geometry GeometryType=\"XYZ\"> \n"
             << "<DataItem Name=\"points\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
             << "field.bin\n"
             << "</DataItem>\n</Geometry>"
             << "<Attribute Name =\"lvlset\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
             << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\"  Endian=\"Little\" Dimensions=\"&Npoints;\">\n"
             << "Phi_t="+ std::to_string(timestep*dt) +".bin\n"
             << "</DataItem></Attribute></Grid>\n";

    delete[] pointCoordinates;
}

double LevelSet::sumLevelSet() {
    double temp = 0;
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++)
		temp = temp + this->at(x, y, z);

    return temp;
}

void LevelSet::initDroplet(array<double, 3> center, double radius, double epsilon) {
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++)
		this->at(x, y, z) = pow(x*dx - center[0], 2) + pow(y*dx - center[1], 2) + pow(z*dx - center[2], 2) - pow(radius, 2);
}

void LevelSet::calculateNextTimestep(double dt) {
    LevelSet tempPhi(*this);

    const array<double, 3> upNormal = {0, 1, 0};
    const array<double, 3> downNormal = {0,-1, 0};
    const array<double, 3> leftNormal = {-1, 0, 0};
    const array<double, 3> rightNormal = {1, 0, 0};
    const array<double, 3> frontNormal = {0, 0, 1};
    const array<double, 3> backNormal = {0, 0, -1};

    for (int x = 0; x < numX; x++) {
    	for (int y = 0; y < numY; y++) {
    		for (int z = 0; z < numZ; z++) {
			//Calculate the flux of the fluid through all sides of the square
			double flux = 0;
			double sp = 0;

			for (int dir = 0; dir < 6; dir++) {
					// no flux across the domain boundary
					if (boundaryCondition == BoundaryCondition::Dirichlet &&
					    ((x == 0 && dir == 2) || (x == numX-1 && dir == 3)
					  || (y == 0 && dir == 1) || (y == numY-1 && dir == 0)
					  || (z == 0 && dir == 5) || (z == numZ-1 && dir == 4)))
							continue;



				switch(dir) {
				case 0:
							sp = field->at((x+1/2)*dx, (y+1)*dx, (z+1/2)*dx) * upNormal;
							if (boundaryCondition == BoundaryCondition::homogeneousNeumann && y == numY - 1) {
								flux += sp*tempPhi.at(x, y, z);
								break;
							}
							flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y+1, z); // upwind flux
							//flux +=(tempPhi.at(x, y+1, z) + tempPhi.at(x, y, z))/2* sp;
							break;
				case 1:
							sp = field->at((x+1/2)*dx, y*dx, (z+1/2)*dx) * downNormal;
							if (boundaryCondition == BoundaryCondition::homogeneousNeumann && y == 0 ) {
								flux += sp*tempPhi.at(x, y, z);
								break;
							}
							flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y-1, z); // upwind flux
							//flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y-1, z))/2 * sp;
							break;
				case 2:
							sp = field->at(x*dx, (y+1/2)*dx, (z+1/2)*dx) * leftNormal;
							if (boundaryCondition == BoundaryCondition::homogeneousNeumann && x == 0) {
								flux += sp*tempPhi.at(x, y, z);
								break;
							}
							flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x-1, y, z); // upwind flux
							//flux += (tempPhi.at(x, y, z) + tempPhi.at(x-1, y, z))/2 * sp;
							break;
				case 3:
							sp = field->at((x+1)*dx, (y+1/2)*dx, (z+1/2)*dx) * rightNormal;
							if (boundaryCondition == BoundaryCondition::homogeneousNeumann && x == numX - 1) {
								flux += sp*tempPhi.at(x, y, z);
								break;
							}
							flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x+1, y, z); // upwind flux
							//flux += (tempPhi.at(x, y, z) + tempPhi.at(x+1, y, z))/2 * sp;
							break;
				case 4:
							if (numZ == 1)
								break;
							sp = field->at((x+1/2)*dx, (y+1/2)*dx, (z+1)*dx) * frontNormal;
							if (boundaryCondition == BoundaryCondition::homogeneousNeumann && z == numZ - 1) {
								flux += sp*tempPhi.at(x, y, z);
								break;
							}
							flux += fmax(sp,0.0)*tempPhi.at(x, y, z+1) + fmin(sp,0.0)*tempPhi.at(x, y, z); // upwind flux
							//flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y, z+1))/2 * sp;
							break;
				case 5:
							if (numZ == 1)
								break;
							sp = field->at((x+1/2)*dx, (y+1/2)*dx, z*dx) * backNormal;
							if (boundaryCondition == BoundaryCondition::homogeneousNeumann && z == 0) {
								flux += sp*tempPhi.at(x, y, z);
								break;
							}
							flux += fmax(sp,0.0)*tempPhi.at(x, y, z-1) + fmin(sp,0.0)*tempPhi.at(x, y, z); // upwind flux
							//flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y, z-1))/2 * sp;
							break;
				}
			}
			this->at(x, y, z) = this->at(x, y, z) - dt/dx*flux;
			}
    	}
    }
}
