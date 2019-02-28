#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include "vecMath2D.hpp"

std::array<double, 2> shearField(double x, double y) {
  return {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y)};
}

void printField(const std::vector <std::vector<double>>& Phi, double epsilon) {
  for (int i = Phi[0].size() - 1; i >= 0; i--) {
    std::ostringstream strs;
    for (int j = 0; j < Phi.size(); j++) {
      if (Phi[j][i] > epsilon) {
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
  std::ofstream outFile("field.csv." + std::to_string(timestep));
  outFile << "X,Y,Z,field" << std::endl;
  for (int y = 0; y < Phi[0].size(); y++)
    for (int x = 0; x < Phi.size(); x++)
      if (abs(Phi[x][y]) < epsilon)
	 outFile << std::to_string(x*dx) + "," + std::to_string(y*dx) +",0," << std::endl;
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

std::vector <std::vector<double>>& initDroplet(std::vector <std::vector<double>>& Phi, std::array<double, 2> center, double radius, double dx, double epsilon) {
  for (int i = 0; i < Phi.size(); i++) {
    for (int j = 0; j < Phi[0].size(); j++) {
      std::array<int, 2> temp = {i, j};
      if (abs(center - temp*dx) < radius + epsilon &&  abs(center - temp*dx) > radius - epsilon) {
	Phi[i][j] = 0.0;
      } else if (abs(center - temp*dx) > radius + epsilon) {
	Phi[i][j] = -1; 
      } else {
	Phi[i][j] = 1;
      }
    }
  }
  return Phi;
}

std::vector <std::vector<double>>& calculateNextTimestep(std::vector <std::vector<double>>& Phi,
							 double dt, double dx,
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
	int integSteps = 100;
	double h = dx/integSteps;
	if ((x == 0 && dir == 2) || (x == Phi.size()-1 && dir == 3) || (y == 0 && dir == 1) || (y == Phi[0].size()-1 && dir == 0))
	  continue;	       
	
	switch(dir) {
	case 0:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = tempPhi[x][y+1]*field(x*dx + (k-1)*h, (y+1)*dx)*upNormal;
	    temp2 = tempPhi[x][y+1]*field(x*dx + k*h, (y+1)*dx)*upNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	case 1:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = tempPhi[x][y-1]*field(x*dx + (k-1)*h, y*dx)*downNormal;
	    temp2 = tempPhi[x][y-1]*field(x*dx + k*h, y*dx)*downNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	case 2:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = tempPhi[x-1][y]*field(x*dx, y*dx + (k-1)*h)*leftNormal;
	    temp2 = tempPhi[x-1][y]*field(x*dx, y*dx +  k*h)*leftNormal;
	    flux += h/2*(temp1 + temp2);
	  }
	  break;
	case 3:
	  for (int k = 1; k <= integSteps; k++) {
	    temp1 = tempPhi[x+1][y]*field((x+1)*dx, y*dx + (k-1)*h)*rightNormal;
	    temp2 = tempPhi[x+1][y]*field((x+1)*dx, y*dx +  k*h)*rightNormal;
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

        int numX, numY, timesteps;
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
	std::array<double, 2> center = {centerX, centerY};
	
	//Initialize empty field
	std::vector< std::vector<double> > Phi(numX);
	for (int i = 0; i < numX; i++) {
	  std::vector<double> temp(numY);
	  Phi[i] = temp;
	}
	
        Phi = initDroplet(Phi, center, radius, dx, 0.005);
	
	double sumAtStart = sumField(Phi);
	for (int i = 0; i < timesteps; i++) {
	  std::cout << "Step " << i << std::endl;
	  //Print field to terminal
	  printField(Phi, 0.3);
	  //Write field to file
	  //writeInterfaceToFile(Phi, dx, 0.3, i);
	  Phi = calculateNextTimestep(Phi, dt, dx, field);
	}
	std::cout << std::endl;
        std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
	std::cout << "Sum of Phi at end: " << sumField(Phi) << std::endl;

       
	return 0;
}

