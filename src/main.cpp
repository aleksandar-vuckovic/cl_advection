/**
 * @file
 *
 * @section Description
 *
 * This file reads all input parameters from the Inputfile and contains the main for-loop which iterates over all timesteps.
 * It opens and writes the output files.
 *
 */


#include <iostream>   // Terminal IO
#include <iomanip>    // IO manipulation, set::fixed, set::setprecision
#include <stdio.h>
#include <string>
#include <vector>
#include "LevelSet.hpp"
#include "VelocityField.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

using std::array;

int main() {

    int numX, numY, numZ;
    double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius, expcpX, expcpY, expcpZ, expAngle, v0, c1, c2, tau, CFL, writestepsFraction;
    numX = numY = numZ = 0;
    lenX = lenY = lenZ = time = centerX = centerY = centerZ = radius = expcpX = expcpY = expcpZ = expAngle = v0 = c1 = c2 = tau = CFL = writestepsFraction = 0;
    bool calculateCurvature = false, writeField = false;
    std::string trackedContactPoint = "left", fieldName = "";
    VelocityField *field = nullptr;

    // Read data from Inputfile
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
		else if (varName == "CFL")
		   CFL = std::stod(value);
		else if (varName == "writestepsFraction")
		    writestepsFraction =std::stod(value);
		else if (varName == "writeField")
			std::stringstream(value) >> std::boolalpha >> writeField;
		else if (varName == "v0")
			v0 = std::stod(value);
		else if (varName == "c1")
			c1 = std::stod(value);
		else if (varName == "c2")
			c2 = std::stod(value);
		else if (varName == "tau")
			tau = std::stod(value);
		else if (varName == "field") {
			fieldName = value;
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
		else if (varName == "trackedContactPoint")
		    trackedContactPoint = value;
		else if (varName == "calculateCurvature")
			std::stringstream(value) >> std::boolalpha >> calculateCurvature;

	    }
	   }
    }

	field = new VelocityField(fieldName, v0, c1, c2, tau, 0, lenX, 0, lenY, 0, lenZ, lenX/numX,lenY/numY,lenZ/numZ);

    double dx = lenX/numX;
    double dy = lenY/numY;
    double dz = lenZ/numZ;

    double dt = 0.0;
    if(lenZ == 1){
      // 2D case
      dt = CFL*std::min(dx,dy)/field->getMaxNormValue();
    }
    else{
      // 3D case
      dt = CFL*std::min(std::min(dx,dy),dz)/field->getMaxNormValue();
    }

    int timesteps = time/dt;
    int writesteps = floor(writestepsFraction*timesteps);
    writesteps = timesteps/(ceil((double)timesteps/writesteps));
    double initCurvature = -1/radius;

    //This is used for the numerical reference solution. Currently this only applies in 2D.
    array<double, 3> n_sigma_init = {-sin(expAngle/180*M_PI), cos(expAngle/180*M_PI), 0};

    LevelSet Phi(numX, numY, numZ, dx, dy, dz, field, trackedContactPoint);

    array<double, 3> center = {centerX, centerY, centerZ};
    array<double, 3> expcp = {expcpX, expcpY, expcpZ};
    std::vector<double> angle(timesteps);
    std::vector<double> curvatureActual(timesteps);
    std::vector<double> curvatureTheoretical(timesteps);

    Phi.initDroplet(center, radius);
    array<double, 3> initCP = Phi.getInitCP(expcp, 0.001);

    int sysRet = system("mkdir data");

    if (sysRet == 256) {
    	std::cout << "Overwriting folder \"data\".\n";
    }
    std::ofstream positionFile("position.csv");
    std::ofstream angleFile("contactAngle.csv");
    std::ofstream curvatureFile("curvature.csv");
    double sumAtStart = Phi.sumLevelSet();

    // XMF file for Paraview
    std::ofstream xmfFile("data/Phi.xmf");

    // main loop

    for (int i = 0; i < timesteps; i++) {
        std::cout << "Step " << i << std::endl;
        //Write field to file
        if (writeField && i % (int)ceil((double)timesteps/writesteps) == 0) {
            Phi.writeToFile(dt, i, timesteps, writesteps, &xmfFile);
        }

		// Calculate the theoretical position of the contact point
        array<double, 3> newCP = Phi.getContactPoint(dt, i, timesteps, initCP);

        // Get the indices of the new contact point
        array<int, 3> newCPIndices = Phi.getContactPointIndices(newCP);

        angle[i] = Phi.getContactAngle(newCPIndices);

        std::cout << "Time: " << i*dt << "\n";
        positionFile << i*dt << ", " << dx*newCPIndices[0] << ", " << newCP[0] << std::endl;


        std::cout << "Actual: " << std::to_string(angle[i]/(2*M_PI)*360) + "\n";
        if (field->getName() == "shearField") {
			double reference_temp = Phi.getReferenceAngleExplicitEuler(dt, i, n_sigma_init, newCP)/M_PI*180;
			angleFile << i*dt << ", " << angle[i]/(2*M_PI)*360 << ", " << reference_temp << "\n";
			std::cout << "Reference: " << reference_temp << "\n";
		} else {
		    double reference_temp = Phi.getReferenceAngleLinearField(i*dt, c1, c2, expAngle/180*M_PI)/M_PI*180.0;
			//double reference_temp = Phi.getReferenceAngleExplicitEuler(dt, i, n_sigma_init, newCP)/M_PI*180;
			angleFile << i*dt << ", " << angle[i]/(2*M_PI)*360 << ", " 	<< reference_temp << "\n";
			std::cout << "Reference: " << reference_temp << "\n";
		}

        if (calculateCurvature) {
            curvatureActual[i] = Phi.getCurvature(newCPIndices);
            curvatureTheoretical[i] = Phi.getReferenceCurvature(dt, i, initCurvature, newCP, newCPIndices);
            curvatureFile << std::to_string(i*dt) + "," + std::to_string(curvatureActual[i]) + "," + std::to_string(curvatureTheoretical[i]) + "\n";

            std::cout << "Measured curvature: " + std::to_string(curvatureActual[i]) + "\n";
            std::cout << "Reference curvature: " + std::to_string(curvatureTheoretical[i])  << std::endl;
		}

        // Calculate numerical flux through all faces of each cell and change Phi accordingly
        Phi.calculateNextTimestep(dt, i);

        positionFile.flush();
		angleFile.flush();
		curvatureFile.flush();
    }
    // Add closing line to the XMF file
    xmfFile << "</Grid>\n</Domain>\n</Xdmf>" << std::endl;

    std::cout << std::endl;
    std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
    std::cout << "Sum of Phi at end: " << Phi.sumLevelSet() << std::endl;

    return 0;
}
