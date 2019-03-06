#include <iostream>
#include <fstream>
#include <sstream>
#include "vecMath3D.hpp"

std::array<double, 3> shearField(double x, double y, double z) {
  return {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y), 0};
}

void writeFieldToFile(const std::vector <std::vector <std::vector<double> >>& Phi, double dx, double epsilon, int timestep) {
  std::ostringstream pointCoordinates;
  std::ostringstream pointPhiValues;
  int Npoints = 0;
  for (int x = 0; x < Phi.size(); x++)
    for (int y = 0; y < Phi[0].size(); y++)
      for (int z = 0; z < Phi[0][0].size(); z++) {
	pointCoordinates << std::to_string(x*dx) + " " + std::to_string(y*dx) + " " + std::to_string(z*dx) + "\n";
	pointPhiValues << std::to_string(Phi[x][y][z]) + "\n"; 
       Npoints++;
     }

  std::ofstream outFile("data/field_t="+ std::to_string(timestep) +".xmf");
  outFile << "<?xml version=\"1.0\" ?>\n"                                                                                         
	  << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"                                                                          
	  << "<!ENTITY Npoints   \" " + std::to_string(Npoints) + "\">\n"                                                    
	  << " ]>"                                                                                                      
          << "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">\n"                                       
          << "<Domain>\n"                                                                                                
	  << "<Grid Name=\"interface\">\n"                                                                            
	  << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"&Npoints;\"/>\n"                      
	  << "<Geometry GeometryType=\"XYZ\">\n"                                                            
	  << "<DataItem Name=\"points\" Format=\"XML\" NumberType=\"Float\" Dimensions=\"&Npoints; 3\">\n";
  
  outFile << pointCoordinates.str();
  outFile <<"</DataItem>\n" <<"</Geometry>\n" << "<Attribute Name=\"Phi\" AttributeType=\"Scalar\" Center=\"Cell\">\n <DataItem Format =\"XML\" NumberType=\"Float\" Dimensions=\"&Npoints;\">\n";
  outFile << pointPhiValues.str();
  outFile << "</DataItem>\n</Attribute>\n</Grid>\n</Domain>\n</Xdmf>\n";
}
  
double sumField(const std::vector <std::vector< std::vector<double> >>& Phi) {
  double temp = 0;
  for (int x = 0; x < Phi.size(); x++) 
    for (int y = 0; y < Phi[0].size(); y++)
      for (int z = 0; z < Phi[0][0].size(); z++)
	temp = temp + Phi[x][y][z];
  
  return temp;
}

std::vector <std::vector<std::vector<double> >>& initDroplet(std::vector <std::vector< std::vector<double> >>& Phi, double dx, std::array<double, 3> center, double radius, double epsilon) {
  for (int x = 0; x < Phi.size(); x++) 
    for (int y = 0; y < Phi[0].size(); y++)
      for (int z = 0; z < Phi[0][0].size(); z++)
	Phi[x][y][z] = pow(x*dx - center[0], 2) + pow(y*dx - center[1], 2) + pow(z*dx - center[2], 2) - pow(radius, 2);

  return Phi;
}

std::vector <std::vector< std::vector<double> >>& calculateNextTimestep(std::vector <std::vector< std::vector<double> >>& Phi,
							 double dx, double dt,
							 std::array<double, 3> (*field) (double x, double y, double z)) {

  std::vector <std::vector< std::vector<double> >> tempPhi(Phi);
  
  const std::array<double, 3> upNormal = {0, 1, 0};
  const std::array<double, 3> downNormal = {0,-1, 0};
  const std::array<double, 3> leftNormal = {-1, 0, 0};
  const std::array<double, 3> rightNormal = {1, 0, 0};
  const std::array<double, 3> frontNormal = {0, 0, 1};
  const std::array<double, 3> backNormal = {0, 0, -1};
  
  for (int x = 0; x < Phi.size(); x++) {
    for (int y = 0; y < Phi[0].size(); y++) {
      for (int z = 0; z < Phi[0][0].size(); z++) {
      //Calculate the flux of the fluid through all sides of the square
	double flux = 0;
      for (int dir = 0; dir < 6; dir++) {
	if ((x == 0 && dir == 2) || (x == Phi.size()-1 && dir == 3)
	    || (y == 0 && dir == 1) || (y == Phi[0].size()-1 && dir == 0)
	    || (z == 0 && dir == 6) || (z = Phi[0][0].size() - 1 && dir == 5))
	  continue;	       
	
	switch(dir) {
	case 0:
	  flux = (tempPhi[x][y][z] + tempPhi[x][y+1][z])/2 * field((x+1/2)*dx, (y+1)*dx, (z+1/2)*dx) * upNormal * dx*dx;
	  break;
	case 1:
	  flux = (tempPhi[x][y][z] + tempPhi[x][y-1][z])/2 * field((x+1/2)*dx, y*dx, (z+1/2)*dx) * downNormal * dx*dx;
	  break;
	case 2:
	  flux = (tempPhi[x][y][z] + tempPhi[x-1][y][z])/2 * field(x*dx, (y+1/2)*dx, (z+1/2)*dx) * leftNormal * dx*dx;
	  break;
	case 3:
	  flux = (tempPhi[x][y][z] + tempPhi[x+1][y][z])/2 * field((x+1)*dx, (y+1/2)*dx, (z+1/2)*dx) * rightNormal * dx*dx;
	  break;
	case 4:
	  flux = (tempPhi[x][y][z] + tempPhi[x][y][z+1])/2 * field((x+1/2)*dx, (y+1/2)*dx, (z+1)*dx) * frontNormal * dx*dx;
	  break;
	case 5:
	  flux = (tempPhi[x][y][z] + tempPhi[x][y][z-1])/2 * field((x+1/2)*dx, (y+1/2)*dx, z*dx) * backNormal * dx*dx;
	  break;
	}
      }
      Phi[x][y][z] = Phi[x][y][z] - dt/(dx*dx)*flux;
      }
    }
  }
  return Phi;
}

int main() {

        int numX, numY, numZ, timesteps, writesteps;
        double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius;
	std::array<double, 3> (*field) (double x, double y, double z);
        
        std::ifstream inFileStream("inputfile");
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
		else if (varName == "field") {
		  if (value == "shearField")
		    field = shearField;
		}
		else if (varName == "centerX")
		  centerX = std::stod(value);
		else if (varName == "centerY")
		  centerY = std::stod(value);
		else if (varName == "centerZ")
		  centerY = std::stod(value);
		else if (varName == "radius")
		  radius = std::stod(value);
	      }
	  }
	}
	    
	

	double dx = lenX/numX;
	double dt = time/timesteps;

	if (dt/dx < 1) {
	  std::cout << "The stability requirement is fullfilled" << std::endl;
	} else {
	  std::cout << "The stability requirement is NOT fullfilled" << std::endl;
	}
	std::array<double, 3> center = {centerX, centerY, centerZ};
	
	//Initialize empty field
	std::vector< std::vector< std::vector<double> >> Phi(numX);
	for (int x = 0; x < numX; x++) {
	  std::vector< std::vector<double> > tempY(numY);
	  for (int y = 0; y < numY; y++) {
	    std::vector<double> tempZ(numZ);
	    tempY[y] = tempZ;
	  }
	  Phi[x] = tempY;
	}
	
        Phi = initDroplet(Phi, dx, center, radius, 0.005);

	system("mkdir data");
	double sumAtStart = sumField(Phi);
	for (int i = 0; i < timesteps; i++) {
	  std::cout << "Step " << i << std::endl;
	  //Write field to file
	  if (i % (timesteps/writesteps) == 0)
	    writeFieldToFile(Phi, dx, 0.01, i);
	  Phi = calculateNextTimestep(Phi, dx, dt, field);
	}
	std::cout << std::endl;
        std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
	std::cout << "Sum of Phi at end: " << sumField(Phi) << std::endl;

       
	return 0;
}

