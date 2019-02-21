#include <iostream>
#include <array>
#include <vector>
#include <ctime>
#include <sstream>

#define _USE_MATH_DEFINES
#include <cmath>

std::array<double, 2> operator+ (std::array<double, 2> vecA, std::array<double, 2> vecB){
  return { vecA[0] + vecB[0], vecA[1] + vecB[1]};
}

std::array<double, 2> operator*(double a, const std::array<double, 2>& vec) {
  return { a*vec[0], a*vec[1]};
}

double operator* (const std::array<double, 2>& vecA, const std::array<double, 2>& vecB) {
  return { vecA[0]*vecB[0] + vecA[1]*vecB[1]};
}


std::array<double, 2> shearField(double x, double y) {
  return {-sin(M_PI*x)*cos(M_PI*y),
	            cos(M_PI*x)*sin(M_PI*y)};
}

void printField(const std::vector <std::vector<double>>& Phi) {
  std::ostringstream strs;
  for (int i = 0; i < Phi.size(); i++) {
    for (int j = 0; j < Phi[0].size(); j++) {
      strs << Phi[i][j];
      strs << "    ";
    }
    strs << std::endl;
  }
  std::cout << strs.str();
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

std::vector <std::vector<double>>& calculateNextTimestep(std::vector <std::vector<double>>& Phi, double dt, double dx) {
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
	    temp1 = Phi[x][y]*shearField(x*dx + (k-1)*h, (y+1)*dx)*upNormal;
	    temp2 = Phi[x][y]*shearField(x*dx + k*h, (y+1)*dx)*upNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  case 1:
	    temp1 = Phi[x][y]*shearField(x*dx + (k-1)*h, y*dx)*downNormal;
	    temp2 = Phi[x][y]*shearField(x*dx + k*h, y*dx)*downNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  case 2:
	    temp1 = Phi[x][y]*shearField(x*dx, y*dx + (k-1)*h)*leftNormal;
	    temp2 = Phi[x][y]*shearField(x*dx, y*dx +  k*h)*leftNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  case 3:
	    temp1 = Phi[x][y]*shearField((x+1)*dx, y*dx + (k-1)*h)*rightNormal;
	    temp2 = Phi[x][y]*shearField((x+1)*dx, y*dx +  k*h)*rightNormal;
	    flux += h/2*(temp1 + temp2);
	    break;
	  }
	}
      }
      Phi[x][y] -= dt/(dx*dx)*flux;
    }
  }
  return Phi;
}

int main() {
	//Number of cells
	int numX = 10;
	int numY = 10;
	double lenX = 1.0;
	double lenY = 1.0;
	double dx = lenX/numX;

	double time = 10.0;
	int timesteps = 1000;
	double dt = time/timesteps;

	auto field = shearField;
	std::vector< std::vector<double> > Phi(numX);

	//Create random field Phi for testing
	srand(std::time(NULL));
	for (int i = 0; i < numX; i++) {
	  std::vector<double> tempY(numY);
	  for (int j = 0; j < numY; j++) {
	    tempY[j] = (rand()%10)/10.0;
	  }
	  Phi[i] = tempY;
	}
	//printField(Phi);
	std::cout << std::endl << sumField(Phi) << std::endl;
	for (int i = 0; i < timesteps; i++)  {
	  Phi = calculateNextTimestep(Phi, dt, dx);
	}
	//printField(Phi);
	std::cout << std::endl;
	std::cout << sumField(Phi) << std::endl;

	return 0;
}

