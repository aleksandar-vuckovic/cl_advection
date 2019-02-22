#include <iostream>
#include <ctime>
#include <sstream>
#include "vecMath2D.hpp"
#include <stdlib.h>

std::array<double, 2> shearField(double x, double y) {
  return {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y)};
}

void printField(const std::vector <std::vector<double>>& Phi, double epsilon) {
  for (int i = Phi[0].size() - 1; i >= 0; i--) {
    std::ostringstream strs;
    for (int j = 0; j < Phi.size(); j++) {
      if (abs(Phi[j][i]) < epsilon) {
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
	Phi[i][j] = -10.0; 
      } else {
	Phi[i][j] = 10.0;
      }
    }
  }
  return Phi;
}

std::vector <std::vector<double>>& calculateNextTimestep(std::vector <std::vector<double>>& Phi, double dt, double dx, auto field) {
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
	for (int k = 1; k <= integSteps; k++) {
	  switch(dir) {
	  case 0:
	    temp1 = Phi[x][y]*field(x*dx + (k-1)*h, (y+1)*dx)*upNormal;
	    temp2 = Phi[x][y]*field(x*dx + k*h, (y+1)*dx)*upNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  case 1:
	    temp1 = Phi[x][y]*field(x*dx + (k-1)*h, y*dx)*downNormal;
	    temp2 = Phi[x][y]*field(x*dx + k*h, y*dx)*downNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  case 2:
	    temp1 = Phi[x][y]*field(x*dx, y*dx + (k-1)*h)*leftNormal;
	    temp2 = Phi[x][y]*field(x*dx, y*dx +  k*h)*leftNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  case 3:
	    temp1 = Phi[x][y]*field((x+1)*dx, y*dx + (k-1)*h)*rightNormal;
	    temp2 = Phi[x][y]*field((x+1)*dx, y*dx +  k*h)*rightNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  }
	}
      }
      Phi[x][y] = Phi[x][y] - dt/(dx*dx)*flux;
    }
  }
  return Phi;
}

int main() {
	//Number of cells
        int numX = 40;
	int numY = 40;
	double lenX = 1.0;
	double lenY = 1.0;
	double dx = lenX/numX;

	//Start time is assumed to be 0s
	double time = 10000;
	int timesteps = 10;
	double dt = time/timesteps;

	auto field = shearField;
	//Initialize empty field
	std::vector< std::vector<double> > Phi(numX);
	for (int i = 0; i < numX; i++) {
	  std::vector<double> temp(numY);
	  Phi[i] = temp;
	}
	std::array<double, 2> center = {0.7, 0.0};
	double radius = 0.2;

	//Create random field Phi for testing
	//	srand(std::time(NULL));
	//for (int i = 0; i < numX; i++) {
	//std::vector<double> tempY(numY);
	//for (int j = 0; j < numY; j++) {
	// tempY[j] = (rand()%100)/100.0;
	//}
	//Phi[i] = tempY;
	//	}

        Phi = initDroplet(Phi, center, radius, dx, 0.02);
	
	std::cout << "Sum of Phi at start: " << sumField(Phi) << std::endl;
	for (int i = 0; i < timesteps; i++) {
	  printField(Phi, 0.0001);
	  Phi = calculateNextTimestep(Phi, dt, dx, field);
	}
	std::cout << std::endl;
	std::cout << "Sum of Phi at end: " << sumField(Phi) << std::endl;

	return 0;
}

