#include <iostream>
#include <fstream>
#include <sstream>
#include "vecMath2D.hpp"

std::array<double, 2> shearField(double x, double y) {
  return {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y)};
}

void printField(const std::vector <std::vector<double>>& Phi, double epsilon) {
  for (int i = Phi[0].size() - 1; i >= 0; i--) {
    std::ostringstream strs;
    for (int j = 0; j < Phi.size(); j++) {
      if (Phi[j][i] < epsilon) {
	strs << "x";
	//strs << " ";
      } else {
	strs << ".";
	//strs << " ";
      }
    }
    std::cout << strs.str();
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void writeInterfaceToFile(const std::vector <std::vector<double>>& Phi, double dx, double epsilon, int timestep) {
  std::ostringstream points;
  int Npoints = 0;
  for (int y = 0; y < Phi[0].size(); y++)
   for (int x = 0; x < Phi.size(); x++)
     if (Phi[x][y] < epsilon) {
       points << std::to_string(x*dx) + " " + std::to_string(y*dx) +" 0.0"                                  << std::endl;  
       Npoints++;
     }

  std::ofstream outFile("data/field_t="+ std::to_string(timestep) +".xmf");
  outFile << "<?xml version=\"1.0\" ?>"                                                                                          << std::endl
	  <<   "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" ["                                                                            << std::endl
	  <<   "<!ENTITY Npoints   \" " + std::to_string(Npoints) + "\">"                                                        << std::endl
	  <<   " ]>"                                                                                                             << std::endl
          <<   "<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">"                                           << std::endl
          <<   "      <Domain>"                                                                                                  << std::endl
	  <<   "          <Grid Name=\"interface\">"                                                                             << std::endl
	  <<   "                  <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"&Npoints;\"/>"                        << std::endl
	  << "                    <Geometry GeometryType=\"XYZ\">"                                                               << std::endl
	  << "                        <DataItem Name=\"points\" Format=\"XML\" NumberType=\"Float\" Dimensions=\"&Npoints; 3\">" << std::endl;
  
  outFile << points.str();
  
  outFile <<"                </DataItem>"                                                                                        << std::endl
          <<"            </Geometry>"                                                                                            << std::endl
	  <<"        </Grid>"                                                                                                    << std::endl
	  <<"    </Domain>"                                                                                                      << std::endl
	  <<"</Xdmf>"                                                                                                            << std::endl;
}
  
double sumField(const std::vector <std::vector<double>>& Phi) {
  double temp = 0;
  for (int i = 0; i < Phi.size(); i++) {
    for (int j = 0; j < Phi[0].size(); j++) {
      temp = temp + Phi[i][j];
    }
  }
  return temp;
}

std::vector <std::vector<double>>& initDroplet(std::vector <std::vector<double>>& Phi, double dx, std::array<double, 2> center, double radius, double epsilon) {
  for (int x = 0; x < Phi.size(); x++) {
    for (int y = 0; y < Phi[0].size(); y++) {
      Phi[x][y] = pow(x*dx-center[0], 2) + pow(y*dx - center[1], 2) - pow(radius, 2);
    }
  }
  return Phi;
}

std::vector <std::vector<double>>& calculateNextTimestep(std::vector <std::vector<double>>& Phi,
							 double dx, double dt,
							 std::array<double, 2> (*field) (double x, double y)) {

  std::vector <std::vector<double>> tempPhi(Phi);
  
  const std::array<double, 2> upNormal = {0,1};
  const std::array<double, 2> downNormal = {0,-1};
  const std::array<double, 2> leftNormal = {-1, 0};
  const std::array<double, 2> rightNormal = {1, 0};
  
  for (int x = 0; x < Phi.size(); x++) {
    for (int y = 0; y < Phi[0].size(); y++) {
      //Calculate the flux of the fluid through all sides of the square
      double flux = 0;
      //Temporary variables for integration
      double temp1;
      double temp2;
      for (int dir = 0; dir < 4; dir++) {
	int integSteps = 2;
	double h = dx/integSteps;
	if ((x == 0 && dir == 2) || (x == Phi.size()-1 && dir == 3) || (y == 0 && dir == 1) || (y == Phi[0].size()-1 && dir == 0))
	  continue;	       
	
	switch(dir) {
	case 0:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = (tempPhi[x][y] + tempPhi[x][y+1])/2*field(x*dx + (k-1)*h, (y+1)*dx)*upNormal;
	    temp2 = (tempPhi[x][y] + tempPhi[x][y+1])/2*field(x*dx + k*h, (y+1)*dx)*upNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	case 1:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = (tempPhi[x][y] + tempPhi[x][y-1])/2*field(x*dx + (k-1)*h, y*dx)*downNormal;
	    temp2 = (tempPhi[x][y] + tempPhi[x][y-1])/2*field(x*dx + k*h, y*dx)*downNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	case 2:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = (tempPhi[x][y] + tempPhi[x-1][y])/2*field(x*dx, y*dx + (k-1)*h)*leftNormal;
	    temp2 = (tempPhi[x][y] + tempPhi[x-1][y])/2*field(x*dx, y*dx +  k*h)*leftNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	case 3:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = (tempPhi[x][y] + tempPhi[x+1][y])/2*field((x+1)*dx, y*dx + (k-1)*h)*rightNormal;
	    temp2 = (tempPhi[x][y] + tempPhi[x+1][y])/2*field((x+1)*dx, y*dx +  k*h)*rightNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	}
	
      }
      Phi[x][y] = Phi[x][y] - dt/(dx*dx)*flux;
    }
  }
  return Phi;
}

int main() {

        int numX, numY, timesteps, writesteps;
        double lenX, lenY, time, centerX, centerY, radius;
	std::array<double, 2> (*field) (double x, double y);
        
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
	        else if (varName == "lenX")
		  lenX = std::stod(value); 
		else if (varName == "lenY")
		  lenY = std::stod(value);
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
	std::array<double, 2> center = {centerX, centerY};
	
	//Initialize empty field
	std::vector< std::vector<double> > Phi(numX);
	for (int i = 0; i < numX; i++) {
	  std::vector<double> temp(numY);
	  Phi[i] = temp;
	}
	
        Phi = initDroplet(Phi, dx, center, radius, 0.005);

	system("mkdir data");
	double sumAtStart = sumField(Phi);
	for (int i = 0; i < timesteps; i++) {
	  std::cout << "Step " << i << std::endl;
	  //Print field to terminal
	  //printField(Phi, 0.01);
	  //Write field to file
	  if (i % (timesteps/writesteps) == 0)
	    writeInterfaceToFile(Phi, dx, 0.01, i);
	  Phi = calculateNextTimestep(Phi, dx, dt, field);
	}
	std::cout << std::endl;
        std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
	std::cout << "Sum of Phi at end: " << sumField(Phi) << std::endl;

       
	return 0;
}

