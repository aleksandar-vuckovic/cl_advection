#include <iostream>
#include <fstream>
#include <sstream>
#include "vecMath3D.hpp"

class LevelSet {
private:
    std::vector<double> data;
    int num_X, num_Y, num_Z;
    double dx;
    std::array<double, 3> (*field) (double x, double y, double z);
    std::array<std::array<double,3>, 3> (*gradientField) (double x, double y, double z);
    
public: 
    LevelSet(int numX, int numY, int numZ, double dx, std::array<double, 3> (*field) (double x, double y, double z),
	     std::array<std::array<double,3>, 3> (*gradientField) (double x, double y, double z)) {
	data = std::vector<double>(numX*numY*numZ);
	num_X = numX;
	num_Y = numY;
	num_Z = numZ;
	this->dx = dx;
	this->field = field;
	this->gradientField = gradientField;
    }
    
    double& at(int x, int y, int z) {
	return data[x + y*num_X + z*num_X*num_Y];
    }
  
    const double& at(int x, int y, int z) const {
	return data[x + y*num_X + z*num_X*num_Y];
    }

    const int& numX() const { return num_X; }
    const int& numY() const { return num_Y; }
    const int& numZ() const { return num_Z; }
    
    std::array<double, 3> getInitCP(double dt, std::array<double, 3> expcp, double epsilon);
    double getContactAngle(double dt, double timestep, int totalTimesteps, std::array<double, 3> initCP);
    double sumLevelSet();
    void writeLevelSetToFile(double epsilon, int timestep);
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
    for (int x = 0; x < this->numX(); x++)
	for (int y = 0; y < this->numY(); y++)
	    for (int z = 0; z < this->numZ(); z++) {
		std::array<double, 3> other = {x*dx, y*dx, z*dx};
		if (abs(candidate - expcp) > abs(other - expcp) && std::abs(this->at(x, y, z)) < epsilon && y == 0) {
		    candidate = other;
		}
	    }
    
    return candidate;
}

double LevelSet::getContactAngle(double dt, double timestep, int totalTimesteps, std::array<double, 3> initCP) {
    std::array<double, 3> &temp = initCP;
    //Calculate current position of contact point
    for (int i = 0; i < timestep; i++)
	temp = temp + dt*field(temp[0], temp[1], temp[2]);

    //Find cell corresponding to this point
    std::array<int, 3> cell = {0, 0, 0};
    for (int x = 0; x < this->numX(); x++)
	for (int y = 0; y < this->numY(); y++)
	    for (int z = 0; z < this->numZ(); z++) {
		std::array<int, 3> other = {x, y, z};
		if (abs(temp - cell*dx) > abs(temp - other*dx))
		    cell = other;
	    }

    //Calculate angle at this cell with finite differences
    double normalX = (this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2]))/(2*dx);
    double normalY = (this->at(cell[0], cell[1] + 1, cell[2]) - this->at(cell[0], cell[1], cell[2]))/dx;
    double normalZ;
    if (this->numZ() > 1)
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
    

void LevelSet::writeLevelSetToFile(double epsilon, int timestep) {
    std::ostringstream pointCoordinates;
    std::ostringstream pointPhiValues;
    int Npoints = 0;
    for (int x = 0; x < this->numX(); x++)
	for (int y = 0; y < this->numY(); y++)
	    for (int z = 0; z < this->numZ(); z++) {
		pointCoordinates << std::to_string(x*dx) + " " + std::to_string(y*dx) + " " + std::to_string(z*dx) + "\n";
		pointPhiValues << std::to_string(this->at(x, y, z)) + "\n"; 
		Npoints++;
	    }

    std::ofstream outFile("data/field_t="+ std::to_string(timestep) +".xmf");
    outFile << "<?xml version=\"1.0\" ?>\n"                                                                                         
	    << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"                                                                          
	    << "<!ENTITY Npoints   \"" + std::to_string(Npoints) + "\">\n"
	    << " ]>"                                                                                                      
	    << "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">\n"                                       
	    << "<Domain>\n"                                                                                                
	    << "<Grid Name=\"interface\">\n"                                                                            
	    << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"&Npoints;\"/>\n"                      
	    << "<Geometry GeometryType=\"XYZ\">\n"                                                            
	    << "<DataItem Name=\"points\" Format=\"XML\" NumberType=\"Float\" Dimensions=\"&Npoints; 3\">\n";
  
    outFile << pointCoordinates.str();
    outFile <<"</DataItem>\n" <<"</Geometry>\n" << "<Attribute Name=\"Phi\" AttributeType=\"Scalar\" Center=\"Cell\">\n<DataItem Format =\"XML\" NumberType=\"Float\" Dimensions=\"&Npoints;\">\n";
    outFile << pointPhiValues.str();
    outFile << "</DataItem>\n</Attribute>\n</Grid>\n</Domain>\n</Xdmf>\n";
}
  
double LevelSet::sumLevelSet() {
    double temp = 0;
    for (int x = 0; x < this->numX(); x++) 
	for (int y = 0; y < this->numY(); y++)
	    for (int z = 0; z < this->numZ(); z++)
		temp = temp + this->at(x, y, z);
  
    return temp;
}

void LevelSet::initDroplet(std::array<double, 3> center, double radius, double epsilon) {
    for (int x = 0; x < this->numX(); x++) 
	for (int y = 0; y < this->numY(); y++)
	    for (int z = 0; z < this->numZ(); z++)
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

    for (int x = 0; x < this->numX(); x++) {
	for (int y = 0; y < this->numY(); y++) {
	    for (int z = 0; z < this->numZ(); z++) {
		//Calculate the flux of the fluid through all sides of the square
		double flux = 0;
		for (int dir = 0; dir < 6; dir++) {
		    if ((x == 0 && dir == 2) || (x == this->numX()-1 && dir == 3)
			|| (y == 0 && dir == 1) || (y == this->numY()-1 && dir == 0)
			|| (z == 0 && dir == 5) || (z == this->numZ()-1 && dir == 4))
			continue;	       
	
		    switch(dir) {
		    case 0:
			flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y+1, z))/2 * field((x+1/2)*dx, (y+1)*dx, (z+1/2)*dx) * upNormal * dx*dx;
			break;
		    case 1:
			flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y-1, z))/2 * field((x+1/2)*dx, y*dx, (z+1/2)*dx) * downNormal * dx*dx;
			break;
		    case 2:
			flux += (tempPhi.at(x, y, z) + tempPhi.at(x-1, y, z))/2 * field(x*dx, (y+1/2)*dx, (z+1/2)*dx) * leftNormal * dx*dx;
			break;
		    case 3:
			flux += (tempPhi.at(x, y, z) + tempPhi.at(x+1, y, z))/2 * field((x+1)*dx, (y+1/2)*dx, (z+1/2)*dx) * rightNormal * dx*dx;
			break;
		    case 4:
			flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y, z+1))/2 * field((x+1/2)*dx, (y+1/2)*dx, (z+1)*dx) * frontNormal * dx*dx;
			break;
		    case 5:
			flux += (tempPhi.at(x, y, z) + tempPhi.at(x, y, z-1))/2 * field((x+1/2)*dx, (y+1/2)*dx, z*dx) * backNormal * dx*dx;
			break;
		    }
		}
		this->at(x, y, z) = this->at(x, y, z) - dt/(dx*dx*dx)*flux;
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
    LevelSet Phi(numX, numY, numZ, dx, field, gradientField);

    if (dt/dx < 1) {
	std::cout << "The stability requirement is fullfilled" << std::endl;
    } else {
	std::cout << "The stability requirement is NOT fullfilled" << std::endl;
    }
    std::array<double, 3> center = {centerX, centerY, centerZ};
    std::array<double, 3> expcp = {expcpX, expcpY, expcpZ};
    std::vector<double> angle(timesteps);
	
    Phi.initDroplet(center, radius, 0.005);
    std::array<double, 3> initCP = Phi.getInitCP(dt, expcp, 0.001);
    
    system("mkdir data");
    std::ofstream angleFile("contactAngleAnalysis.txt");
    double sumAtStart = Phi.sumLevelSet();
    for (int i = 0; i < timesteps; i++) {
	std::cout << "Step " << i << std::endl;
	//Write field to file
	if (i % (timesteps/writesteps) == 0) {
	    Phi.writeLevelSetToFile(0.01, i);
	}
	angle[i] = Phi.getContactAngle(dt, i, timesteps, initCP);
	angleFile << std::to_string(i*dt) + " " + std::to_string(angle[i]/(2*M_PI)*360) + "\n";
	std::cout  << std::to_string(i*dt) + " " + std::to_string(angle[i]/(2*M_PI)*360) + "\n";
	Phi.calculateNextTimestep(dt);
    }
    std::cout << std::endl;
    std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
    std::cout << "Sum of Phi at end: " << Phi.sumLevelSet() << std::endl;
    
    return 0;
}

