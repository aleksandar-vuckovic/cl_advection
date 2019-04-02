#include <iostream>   // Terminal IO
#include <iomanip>    // IO manipulation, set::fixed, set::setprecision
#include <fstream>    // Filestream
#include <sstream>    // Stringstream
#include <stdio.h>
#include <string>
#include "vecMath3D.hpp"

template <class T>
class Field {
private:
    std::vector<T> data;

protected:
    int numX, numY, numZ;

public:
    Field(int numX, int numY, int numZ) {
        data = std::vector<T>(numX*numY*numZ);
        this->numX = numX;
        this->numY = numY;
        this->numZ = numZ;
    }

    T& at(int x, int y, int z) {
	return data[x + y*numX + z*numX*numY];
    }

    const T& at(int x, int y, int z) const {
	return data[x + y*numX + z*numX*numY];
    }
};

class LevelSet : Field<double> {
private:
    double dx;
    std::array<double, 3> (*field) (double x, double y, double z);
    std::array<std::array<double,3>, 3> (*gradientField) (double x, double y, double z);

public:
    LevelSet(int numX, int numY, int numZ, double dx, std::array<double, 3> (*field) (double x, double y, double z),
	     std::array<std::array<double,3>, 3> (*gradientField) (double x, double y, double z)) : Field<double>(numX, numY, numZ) {
	this->dx = dx;
	this->field = field;
	this->gradientField = gradientField;
    }

    std::array<double, 3> getInitCP(double dt, std::array<double, 3> expcp, double epsilon);
    double getContactAngle(double dt, double timestep, std::array<double, 3> initCP);
    double getReferenceCurvature(double dt, double timestep, double initCurvature, std::array<double, 3> initCP) const;
    //For now, getCurvature() only works for the stationary droplet
    double getCurvature(double dt, int timestep,  std::array<double, 3> init) const;
    double sumLevelSet();
    void writeToFile(double epsilon, double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile);
    void initDroplet(std::array<double, 3> center, double radius, double epsilon);
    void calculateNextTimestep(double dt);
};

std::array<double, 3> shearField(double x, double y, double z) {
    return {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y), 0};
}

std::array<std::array<double, 3>, 3> gradShearField(double x, double y, double z) {
    std::array<std::array<double, 3>, 3> tempReturn;
    tempReturn[0] = {-M_PI*cos(M_PI*x)*cos(M_PI*y), -M_PI*sin(M_PI*x)*sin(M_PI*y), 0};
    tempReturn[1] = {-M_PI*sin(M_PI*x)*sin(M_PI*y), M_PI*cos(M_PI*x)*cos(M_PI*y), 0};
    tempReturn[2] = {0, 0, 0};
    return tempReturn;
}

std::array<double, 3> LevelSet::getInitCP(double dt, std::array<double, 3> expcp, double epsilon) {
    std::array<double, 3> candidate = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++) {
		std::array<double, 3> other = {x*dx, y*dx, z*dx};
		if (abs(candidate - expcp) > abs(other - expcp) && std::abs(this->at(x, y, z)) < epsilon && y == 0) {
		    candidate = other;
		}
	    }

    return candidate;
}

double LevelSet::getContactAngle(double dt, double timestep, std::array<double, 3> initCP) {
    std::array<double, 3> &temp = initCP;
    //Calculate current position of contact point
    for (int i = 0; i < timestep; i++)
	temp = temp + dt*field(temp[0], temp[1], temp[2]);

    //Find cell corresponding to this point
    std::array<int, 3> cell = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++) {
		std::array<int, 3> other = {x, y, z};
		if (abs(temp - cell*dx) > abs(temp - other*dx))
		    cell = other;
	    }

    //Calculate angle at this cell with finite differences
    double normalX = (this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2]))/(2*dx);
    //double normalY = (this->at(cell[0], cell[1] + 1, cell[2]) - this->at(cell[0], cell[1], cell[2]))/dx; // first order difference quotient
    double normalY = (-this->at(cell[0], cell[1] + 2, cell[2]) + 4.0*this->at(cell[0], cell[1] + 1, cell[2]) - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dx); // second order difference quotient
    double normalZ;
    if (this->numZ > 1)
	normalZ = (this->at(cell[0], cell[1], cell[2]+1) - this->at(cell[0]-1, cell[1], cell[2]-1))/(2*dx);
    else
	normalZ = 0;
    std::array<double ,3> normal = {normalX, normalY, normalZ};
    normal = normal/abs(normal);
    return acos(normal[1]);
}

/* double getReferenceAngle(double dt, double timestep, int totalTimesteps, double contactAngle) {
    //Calculate current position of contact point
    for (int i = 0; i < timestep; i++)
	temp = temp + dt*field(temp[0], temp[1], temp[2]);
*/

/*
  Currently this only works for the shear field case with theta == pi/2.
 */
double LevelSet::getReferenceCurvature(double dt, double timestep, double initCurvature, std::array<double, 3> initCP) const {
    std::array<double, 3> &CP = initCP;
    //Calculate current position of contact point
    for (int i = 0; i < timestep; i++)
	CP = CP + dt*field(CP[0], CP[1], CP[2]);

     /* Find cell corresponding to this point
       The array "cell" contains the coordinates of the cell */
    std::array<int, 3> cell = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++) {
		std::array<int, 3> other = {x, y, z};
		if (abs(CP - cell*dx) > abs(CP - other*dx))
		    cell = other;
	    }
    
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
        std::array<double, 3> normal = {normalX, normalY, normalZ};
        normal = normal/abs(normal);

        // This is the second derivative of v in the tau direction (= y direction)
        std::array<double, 3> temp = {M_PI*M_PI*sin(M_PI*CP[0])*cos(M_PI*CP[1]), -M_PI*M_PI*cos(M_PI*CP[0])*sin(M_PI*CP[1]), 0};
        
        std::array<double,3> tau = {normal[1], -normal[0], 0};
        
        curvature = curvature + dt*(temp*normal - 2*curvature*(gradientField(CP[0], CP[1], CP[2])*tau)*tau);
    }

    return curvature;
}

double LevelSet::getCurvature(double dt, int timestep, std::array<double, 3> init) const {
    //Calculate current position of contact point
    for (int i = 0; i < timestep; i++)
	init = init + dt*field(init[0], init[1], init[2]);
    
    /* Find cell corresponding to this point
       The array "cell" contains the coordinates of the cell */
    std::array<int, 3> cell = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++) {
		std::array<int, 3> other = {x, y, z};
		if (abs(init - cell*dx) > abs(init - other*dx))
		    cell = other;
	    }

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
    Field<std::array<double, 3> > localField(sidelength, sidelength, sidelengthZ);

    for (int x = -sidelength/2; x <= sidelength/2; x++)
        for (int y = 0; y <= sidelength/2; y++)
            for (int z = -sidelength/2; z <= sidelength/2; z++) {
                if (this->numZ == 1 && z != 0) {
                    continue;
                }
                std::array<int, 3> temp = {x, y, z};
                temp = temp + cell;

                double normalX = (this->at(temp[0] + 1, temp[1], temp[2]) - this->at(temp[0]-1, temp[1], temp[2])) / (2*dx);
                // second order difference quotient
                double normalY = (-this->at(temp[0], temp[1] + 2, temp[2]) + 4.0*this->at(temp[0], temp[1] + 1, temp[2]) - 3.0*this->at(temp[0], temp[1], temp[2])) / (2*dx);
                double normalZ;
                if (this->numZ > 1)
                    normalZ = (this->at(temp[0], temp[1], temp[2] + 1) - this->at(temp[0], temp[1], temp[2] - 1)) / (2*dx);
                else
                    normalZ = 0;
                std::array<double ,3> normal = {normalX, normalY, normalZ};
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
    std::array<double, 3> normal = localField.at(sidelength/2, sidelength/2, sidelengthZ/2);
    std::array<double, 3> tau = {normal[1], -normal[0], 0};
    std::array<double, 3> row1 = {dnx_dx, dnx_dy, 0};
    std::array<double, 3> row2 = {dny_dx, dny_dy, 0};
    std::array<double, 3> row3 = {0,      0,      0};
    std::array< std::array<double, 3>, 3> gradNormal = {row1, row2, row3};

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
                    pointCoordinates[index] = y*dx;
                    pointCoordinates[index +1] = x*dx;
                    pointCoordinates[index +2] = z*dx;
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

void LevelSet::initDroplet(std::array<double, 3> center, double radius, double epsilon) {
    for (int x = 0; x < this->numX; x++)
	for (int y = 0; y < this->numY; y++)
	    for (int z = 0; z < this->numZ; z++)
		this->at(x, y, z) = pow(x*dx - center[0], 2) + pow(y*dx - center[1], 2) + pow(z*dx - center[2], 2) - pow(radius, 2);
}

void LevelSet::calculateNextTimestep(double dt) {

    LevelSet tempPhi(*this);

    const std::array<double, 3> upNormal = {0, 1, 0};
    const std::array<double, 3> downNormal = {0,-1, 0};
    const std::array<double, 3> leftNormal = {-1, 0, 0};
    const std::array<double, 3> rightNormal = {1, 0, 0};
    const std::array<double, 3> frontNormal = {0, 0, 1};
    const std::array<double, 3> backNormal = {0, 0, -1};

    for (int x = 0; x < this->numX; x++) {
	for (int y = 0; y < this->numY; y++) {
	    for (int z = 0; z < this->numZ; z++) {
		//Calculate the flux of the fluid through all sides of the square
		double flux = 0;
                double sp = 0;

		for (int dir = 0; dir < 6; dir++) {
		    if ((x == 0 && dir == 2) || (x == this->numX-1 && dir == 3)
			|| (y == 0 && dir == 1) || (y == this->numY-1 && dir == 0)
			|| (z == 0 && dir == 5) || (z == this->numZ-1 && dir == 4))
			continue; // no flux across the domain boundary

		    switch(dir) {
		    case 0:
                        sp = field((x+1/2)*dx, (y+1)*dx, (z+1/2)*dx) * upNormal;
                        flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y+1, z); // upwind flux
                        //flux +=(tempPhi.at(x, y+1, z) + tempPhi.at(x, y, z))/2* sp;
                        break;
 		    case 1:
                        sp = field((x+1/2)*dx, y*dx, (z+1/2)*dx) * downNormal;
                        flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y-1, z); // upwind flux
			//flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y-1, z))/2 * sp;
			break;
		    case 2:
                        sp = field(x*dx, (y+1/2)*dx, (z+1/2)*dx) * leftNormal;
                        flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x-1, y, z); // upwind flux
			//flux += (tempPhi.at(x, y, z) + tempPhi.at(x-1, y, z))/2 * sp;
			break;
		    case 3:
                        sp = field((x+1)*dx, (y+1/2)*dx, (z+1/2)*dx) * rightNormal;
                        flux += fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x+1, y, z); // upwind flux
			//flux += (tempPhi.at(x, y, z) + tempPhi.at(x+1, y, z))/2 * sp;
			break;
		    case 4:
                        sp = field((x+1/2)*dx, (y+1/2)*dx, (z+1)*dx) * frontNormal;
                        flux += fmax(sp,0.0)*tempPhi.at(x, y, z+1) + fmin(sp,0.0)*tempPhi.at(x, y, z); // upwind flux
			//flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y, z+1))/2 * sp;
			break;
		    case 5:
                        sp = field((x+1/2)*dx, (y+1/2)*dx, z*dx) * backNormal;
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

int main() {

    int numX, numY, numZ, timesteps, writesteps, numCores;
    double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius, expcpX, expcpY, expcpZ, expAngle;
    std::array<double, 3> (*field) (double x, double y, double z);
    std::array<std::array<double, 3>, 3> (*gradientField) (double x, double y, double z);

    std::ifstream inFileStream("Inputfile");
    std::string line, varName, value;

    while(std::getline(inFileStream, line)) {
	std::istringstream linestream(line);
	if(std::getline(linestream, varName, '=')) {
	    if ( std::getline(linestream, value)) {
	        if (varName == "numX")
		    numX = std::stoi(value);
		else if (varName == "numY")
		    numY = std::stoi(value);
		else if (varName == "numZ")
		    numZ = std::stoi(value);
	        else if (varName == "lenX")
		    lenX = std::stod(value);
		else if (varName == "lenY")
		    lenY = std::stod(value);
		else if (varName == "lenZ")
		    lenZ = std::stod(value);
		else if (varName == "time")
		    time = std::stod(value);
		else if (varName == "timesteps")
		    timesteps = std::stoi(value);
		else if (varName == "writesteps")
		    writesteps =std::stoi(value);
		else if (varName == "numCores")
		    numCores = std::stoi(value);
		else if (varName == "field") {
		    if (value == "shearField") {
			field = shearField;
			gradientField = gradShearField;
		    }
		}
		else if (varName == "centerX")
		    centerX = std::stod(value);
		else if (varName == "centerY")
		    centerY = std::stod(value);
		else if (varName == "centerZ")
		    centerZ = std::stod(value);
		else if (varName == "radius")
		    radius = std::stod(value);
		else if (varName == "expcpX")
		    expcpX = std::stod(value);
		else if (varName == "expcpY")
		    expcpY = std::stod(value);
		else if (varName == "expcpZ")
		    expcpZ = std::stod(value);
		else if (varName == "expAngle")
		    expAngle = std::stod(value);

	    }
	}
    }

    double dx = lenX/numX;
    double dt = time/timesteps;
    double initCurvature = -1/radius;
    LevelSet Phi(numX, numY, numZ, dx, field, gradientField);

    if (dt/dx < 1) {
	std::cout << "The stability requirement is fullfilled" << std::endl;
    } else {
	std::cout << "The stability requirement is NOT fullfilled" << std::endl;
    }
    std::array<double, 3> center = {centerX, centerY, centerZ};
    std::array<double, 3> expcp = {expcpX, expcpY, expcpZ};
    std::vector<double> angle(timesteps);
    std::vector<double> curvatureActual(timesteps);
    std::vector<double> curvatureTheoretical(timesteps);

    Phi.initDroplet(center, radius, 0.005);
    std::array<double, 3> initCP = Phi.getInitCP(dt, expcp, 0.001);

    system("mkdir data");
    std::ofstream angleFile("contactAngle.csv");
    std::ofstream curvatureFile("curvature.csv");
    double sumAtStart = Phi.sumLevelSet();

    // XMF file for Paraview
    std::ofstream xmfFile("data/Phi.xmf");

    for (int i = 0; i < timesteps; i++) {
	std::cout << "Step " << i << std::endl;
	//Write field to file
	if (i % (timesteps/writesteps) == 0) {
	    Phi.writeToFile(0.01, dt, i, timesteps, writesteps, &xmfFile);
	}
	angle[i] = Phi.getContactAngle(dt, i, initCP);
        std::cout << "Time: " + std::to_string(i*dt) + "\n"; 
	angleFile << std::to_string(i*dt) + ", " + std::to_string(angle[i]/(2*M_PI)*360) + "\n";
	std::cout << std::to_string(angle[i]/(2*M_PI)*360) + "\n";

        curvatureActual[i] = Phi.getCurvature(dt, i, initCP);
        curvatureTheoretical[i] = Phi.getReferenceCurvature(dt, i, initCurvature, initCP);
        curvatureFile << std::to_string(i*dt) + "," + std::to_string(curvatureActual[i]) + "," + std::to_string(curvatureTheoretical[i]) + "\n";
        std::cout << "Measured curvature: " + std::to_string(curvatureActual[i]) + "\n";
        std::cout << "Reference curvature:" + std::to_string(curvatureTheoretical[i])  << std::endl;
        
	Phi.calculateNextTimestep(dt);
    }
    xmfFile << "</Grid>\n</Domain>\n</Xdmf>" << std::endl;
    std::cout << std::endl;
    std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
    std::cout << "Sum of Phi at end: " << Phi.sumLevelSet() << std::endl;

    return 0;
}

