#include <iostream>
#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>

class unit{
	public:
		//The neighboring volumes
		unit right, left, top, bottom;

		//The position in the plane
		double posX;
		double posY;
};

std::vector<double> shearField(double x, double y){
	std::vector<double> temp;
	temp.push_back(-sin(M_PI*x)*cos(M_PI*y));
	temp.push_back(cos(M_PI*x)*sin(M_PI*y));
	return temp;
}

void calculateNextTimestep(double Phi, double dt, double dx, int numX, int numY) {
	const std::vector<double> upNormal = {0,1};
	const std::vector<double> downNormal = {0,-1};
	const std::vector<double> rightNormal = {1, 0};
	const std::vector<double> leftNormal = {-1, 0};

	for (int x = 0; x < numX; x++){
		for (int y = 0; y < numY; y++){
			//Iterate over all 4 sides of the square
			for (int dir = 0; dir < 4; dir++){
				//Calculate the integral
				double integ = 0;
				int integSteps = 100;
				double h = dx/integSteps;
				for (int i = 0; i < integSteps; i++){

				}


			}
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

	double Phi[numX][numY];



	return 0;
}
