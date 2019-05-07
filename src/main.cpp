#include <iostream>   // Terminal IO
#include <iomanip>    // IO manipulation, set::fixed, set::setprecision
#include <stdio.h>
#include <string>
#include <vector>
#include "LevelSet.hpp"
#include "VelocityField.hpp"
#include "BoundaryCondition.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

using std::array;

int main() {

    int numX, numY, numZ, timesteps, writesteps, numCores;
    double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius, expcpX, expcpY, expcpZ, expAngle, v0, c1, c2;
    bool calculateCurvature, writeVOF;
    VelocityField *field;
    BoundaryCondition boundaryCondition;

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
		else if (varName == "writeVOF")
			std::stringstream(value) >> std::boolalpha >> writeVOF;
		else if (varName == "numCores")
		    numCores = std::stoi(value);
		else if (varName == "v0") {
			v0 = std::stod(value);
		}
		else if (varName == "c1") {
			c1 = std::stod(value);
		}
		else if (varName == "c2") {
			c2 = std::stod(value);
		}
		else if (varName == "field") {
			field = new VelocityField(value, v0, c1, c2);
		}
		else if (varName == "BoundaryCondition") {
			if (value == "homogeneousNeumann")
				boundaryCondition = BoundaryCondition::homogeneousNeumann;
			else if (value == "Dirichlet")
				boundaryCondition = BoundaryCondition::Dirichlet;
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
		else if (varName == "calculateCurvature")
			std::stringstream(value) >> std::boolalpha >> calculateCurvature;

	    }
	}
    }

    double dx = lenX/numX;
    double dt = time/timesteps;
    double initCurvature = -1/radius;
    LevelSet Phi(numX, numY, numZ, dx, field, boundaryCondition);

    if (dt/dx < 1) { // TODO this criterion does only make sense for |v|=1
	std::cout << "The stability requirement is fullfilled" << std::endl;
    } else {
	std::cout << "The stability requirement is NOT fullfilled" << std::endl;
    }
    array<double, 3> center = {centerX, centerY, centerZ};
    array<double, 3> expcp = {expcpX, expcpY, expcpZ};
    std::vector<double> angle(timesteps);
    std::vector<double> curvatureActual(timesteps);
    std::vector<double> curvatureTheoretical(timesteps);

    Phi.initDroplet(center, radius, 0.005);
    array<double, 3> initCP = Phi.getInitCP(dt, expcp, 0.001);

    system("mkdir data");
    std::ofstream positionFile("position.csv");
    std::ofstream angleFile("contactAngle.csv");
    std::ofstream curvatureFile("curvature.csv");
    double sumAtStart = Phi.sumLevelSet();

    // XMF file for Paraview
    std::ofstream xmfFile("data/Phi.xmf");

    for (int i = 0; i < timesteps; i++) {
		std::cout << "Step " << i << std::endl;
		//Write field to file
		if (writeVOF && i % (timesteps/writesteps) == 0) {
			Phi.writeToFile(0.01, dt, i, timesteps, writesteps, &xmfFile);
		}
		array<double, 3> newCP = Phi.getContactPoint(dt, i, timesteps, initCP);
		array<int, 3> newCPCoord = Phi.getContactPointCoordinates(newCP);
    angle[i] = Phi.getContactAngle(dt, i, newCPCoord);
		std::cout << "Time: " + std::to_string(i*dt) + "\n";
		positionFile << std::to_string(i*dt) + ", " << std::to_string(dx*newCPCoord[0])  + ", "<< std::to_string(newCP[0]) << std::endl;
		angleFile << std::to_string(i*dt) + ", " + std::to_string(angle[i]/(2*M_PI)*360)
                    + std::to_string(Phi.getReferenceAngleLinearField(i*dt, c1, c2, expAngle/360*(2*M_PI))/(2*M_PI)*360) + "\n";
		std::cout << "Actual: " << std::to_string(angle[i]/(2*M_PI)*360) + "\n";
                std::cout << "Reference: " << std::to_string(Phi.getReferenceAngleLinearField(i*dt, c1, c2, expAngle/360*(2*M_PI))/(2*M_PI)*360) + "\n";

		if (calculateCurvature) {
			curvatureActual[i] = Phi.getCurvature(dt, i, newCPCoord);
			curvatureTheoretical[i] = Phi.getReferenceCurvature(dt, i, initCurvature, newCP, newCPCoord);
			curvatureFile << std::to_string(i*dt) + "," + std::to_string(curvatureActual[i]) + "," + std::to_string(curvatureTheoretical[i]) + "\n";

			std::cout << "Measured curvature: " + std::to_string(curvatureActual[i]) + "\n";
			std::cout << "Reference curvature: " + std::to_string(curvatureTheoretical[i])  << std::endl;
		}
                // Calculate numerical flux through all faces of each cell and change Phi accordingly
		Phi.calculateNextTimestep(dt);

		positionFile.flush();
		angleFile.flush();
		curvatureFile.flush();
    }
    xmfFile << "</Grid>\n</Domain>\n</Xdmf>" << std::endl;
    std::cout << std::endl;
    std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
    std::cout << "Sum of Phi at end: " << Phi.sumLevelSet() << std::endl;

    return 0;
}
