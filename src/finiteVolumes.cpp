#include <iostream>
#include <array>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <cmath>

std::array<double, 2> operator*(double a, std::array<double, 2> vec) {
  return { vec[0]*a, vec[1]*a};
}

double operator* (std::array<double, 2> vecA, std::array<double, 2> vecB) {
  return { vecA[0]*vecB[0] + vecA[1]*vecB[1]};
}


std::array<double, 2> shearField(double x, double y){
	return {-sin(M_PI*x)*cos(M_PI*y),
			 cos(M_PI*x)*sin(M_PI*y)};
}

void calculateNextTimestep(std::vector <std::vector<double>> Phi, double dt, double dx, int numX, int numY) {
	const std::array<double, 2> upNormal = {0,1};
	const std::array<double, 2> downNormal = {0,-1};
	const std::array<double, 2> rightNormal = {1, 0};
	const std::array<double, 2> leftNormal = {-1, 0};

	for (int x = 0; x < numX; x++){
		for (int y = 0; y < numY; y++){
			//Calculate the flux of the fluid through all sides of the square
			double flux = 0;
			//Temporary variables for integration
			std::array<double, 2> temp1;
			std::array<double, 2> temp2;
			for (int dir = 0; dir < 4; dir++){
				int integSteps = 100;
				double h = dx/integSteps;
				for (int k = 1; k < integSteps; k++){
					switch(dir) {
						case 1:
							temp1 = Phi[x][y]*shearField(x*dx + (k-1)*h, y*dx);
							temp2 = Phi[x][y]*shearField(x*dx + k*h, y*dx);
							flux += h/2*temp1*temp2;
							break;
						case 2:
							temp1 = Phi[x][y]*shearField((x+1)*dx + (k-1)*h, (y+1)*dx);
							temp2 = Phi[x][y]*shearField(x*dx + k*h, (y+1)*dx);
							flux += h/2*temp1*temp2;
							break;
						case 3:
							temp1 = Phi[x][y]*shearField(x*dx, y*dx + (k-1)*h);
							temp2 = Phi[x][y]*shearField(x*dx, y*dx +  k*h);
							flux += h/2*temp1*temp2;
							break;
						case 4:
							temp1 = Phi[x][y]*shearField((x+1)*dx, y*dx + (k-1)*h);
							temp2 = Phi[x][y]*shearField((x+1)*dx, y*dx +  k*h);
							flux += h/2*temp1*temp2;
							break;

					}
				}
			}
			Phi[x][y] -= dt/(dx*dx)*flux;
		}
	}
}

int main() {
	//Number of cells
	int numX = 100;
	int numY = 100;
	double lenX = 1.0;
	double lenY = 1.0;
	double dx = lenX/numX;

	double time = 1.0;
	int timesteps = 1000;
	double dt = time/timesteps;

	auto field = shearField;
+-
	std::vector< std::vector<double> > Phi;

	//test simulation


	std::cout << rand();


	return 0;
}
