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
#include <chrono>
#include "LevelSet.hpp"
#include "VelocityField.hpp"
#include "Streamlines.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

using std::array;

/**
 * The main function of the program.
 *
 * Reads the Inputfile and handles the top-level loop which evolves the field with time.
 * The Inputfile parameters are as follows:
 * numX, numY, numZ | Number of cells in the given direction
 * lenX, lenY, lenZ | Length of simulation plane in the given direction
 * time | The simulation time
 * CFL  | Courant-Friedrics-Lewy-number
 * writestepsFraction | Fraction of total timesteps that will be written to disk
 * writeField | Whether to write the levelset field binary files to disk. The files contactAngle.csv and position.csv are always written (see below)
 * v0, c1, c2, field |  The field and its parameters. For the navier field all three parameters are relevant, but for the shear field, c1 and c2 are ignored and only the value of v0 matters.
 * centerX, centerY, centerZ | The initial center of the droplet.
 * expcpX, expcpY, expcpZ, expAngle | The expected coordinates for the initial contact point and the expected initial contact angle. Those are used for the reference solutions.
 * trackedContactPoint | Whether to track the left or right contact point. Only applicable in 2D.
 */
int main() {
	auto start = std::chrono::system_clock::now();

    int numX, numY, numZ, threads;
    double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius, expcpX, expcpY, expcpZ, expAngle, v0, c1, c2, c3, tau, CFL, writestepsFraction;
    numX = numY = numZ = 0;
    threads = 1;
    lenX = lenY = lenZ = time = centerX = centerY = centerZ = radius = expcpX = expcpY = expcpZ = expAngle = v0 = c1 = c2 = c3 = tau = CFL = writestepsFraction = 0;
    bool writeField = false, calculateCurvature = false;
    std::string trackedContactPoint = "left", fieldName = "";
    VelocityField *field = nullptr;

    // Read data from Inputfile
    std::ifstream inFileStream("Inputfile");
    if(!inFileStream.good()){
        std::cout << "Error: File 'Inputfile' not found in current folder.\n";
        exit(-1);
    }

    std::string line, varName, value;

  while(std::getline(inFileStream, line)) {
	std::istringstream linestream(line);
	if(std::getline(linestream, varName, '=')) {
	    if (std::getline(linestream, value)) {
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
		else if (varName == "c3")
		    c3 = std::stod(value);
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
		else if (varName == "threads")
			threads = std::stoi(value);
	    }
	   }
    }

  	//Parallel computing
  	omp_set_num_threads(threads);

    field = new VelocityField(fieldName, v0, c1, c2, c3, tau, 0, lenX, 0, lenY, 0, lenZ, lenX/numX,lenY/numY,lenZ/numZ);

    double dx = lenX/numX;
    double dy = lenY/numY;
    double dz = lenZ/numZ;
    double initCurvature = -1/radius;

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

    //This is used for the numerical reference solution. Currently this only applies in 2D.
    array<double, 3> n_sigma_init = {-sin(expAngle/180*M_PI), cos(expAngle/180*M_PI), 0};


    array<double, 3> center = {centerX, centerY, centerZ};
    array<double, 3> expcp = {expcpX, expcpY, expcpZ};
    std::vector< array<double, 3> > positionTheoretical(timesteps);
    std::vector<double> angle(timesteps);
    std::vector<double> curvatureActualDivergence(timesteps);
    std::vector<double> curvatureActualHeight(timesteps);
    std::vector<double> curvatureTheoretical(timesteps);


    LevelSet Phi(numX, numY, numZ, dx, dy, dz, field, trackedContactPoint, &positionTheoretical);
    Streamlines streamlines(numX, numY, numZ, *field, dt);

    Phi.initDroplet(center, radius);
    array<double, 3> initCP = Phi.getInitCP(expcp, 0.001);

    int sysRet = system("mkdir data");

    if (sysRet == 256) {
    	std::cout << "Overwriting folder \"data\".\n";
    }
    std::ofstream positionFile("position.csv");
    std::ofstream angleFile("contactAngle.csv");
    std::ofstream curvatureFile;
    if (calculateCurvature)
        curvatureFile = std::ofstream("curvature.csv");
    double sumAtStart = Phi.sumLevelSet();

    // XMF file for Paraview
    std::ofstream xmfFile("data/Phi.xmf");

    // Loop to calculate reference data of position
    for (int i = 0; i < timesteps; ++i) {
        positionTheoretical[i] = Phi.getContactPointExplicitEuler(dt, i, initCP);
    }

    // main loop
    for (int i = 0; i < timesteps; ++i) {
        std::cout << "Step " << i << std::endl;
        //Write field to file
        if (writeField && i % (int)ceil((double)timesteps/writesteps) == 0) {
            if (i == 0)
                streamlines.writeToFile();
            Phi.writeToFile(dt, i, timesteps, writesteps, &xmfFile);
        }

	// Calculate the reference position of the contact point
        array<double, 3> newCPReference;
        if (field->getName() == "navierField" || field->getName() == "timeDependentNavierField")
        	newCPReference = Phi.getContactPointLinearField(dt*i, c1, expcpX, v0);
        else
            newCPReference = positionTheoretical[i];

        // Get the indices of the new contact point
        array<int, 3> newCPIndicesActual = Phi.getContactPointIndices(newCPReference);
        array<double, 3> newCPActual = Phi.getContactPoint(newCPIndicesActual);

        // Evaluate te Contact Angle numerically based on Phi
        angle[i] = Phi.getContactAngle(newCPIndicesActual);

        // Output to command line and positionFile
        std::cout << "Time: " << std::to_string(i*dt) << "\n";
        positionFile << std::to_string(i*dt) << ", "
                     << std::to_string(newCPActual[0]) << ", "
                     << std::to_string(newCPReference[0]) << std::endl;


        std::cout << "Actual: " << std::to_string(angle[i]/(2*M_PI)*360) + "\n";
        if (field->getName() == "navierField" || field->getName() == "timeDependentLinearField") {
            double reference_temp = Phi.getReferenceAngleLinearField(i*dt, c1, c2, expAngle/180*M_PI)/M_PI*180;
			angleFile << std::to_string(i*dt) << ", " << angle[i]/(2*M_PI)*360 << ", " 	<< reference_temp << "\n";
			std::cout << "Reference: " << reference_temp << "\n";
		} else {
            // compute reference for the contact angle (numerically)
            double reference_temp = Phi.getReferenceAngleExplicitEuler(dt, i, n_sigma_init, initCP)/M_PI*180;
            angleFile << std::to_string(i*dt) << ", " << angle[i]/(2*M_PI)*360 << ", " << reference_temp << "\n";
            std::cout << "Reference: " << reference_temp << "\n";
		}
        
        if (calculateCurvature) {
            curvatureActualDivergence[i] = Phi.getCurvatureDivergence(newCPIndicesActual);
            curvatureActualHeight[i] = Phi.getCurvatureHeight(newCPIndicesActual);

            if (field->getName() == "quadraticField") {
                curvatureTheoretical[i] = Phi.getReferenceCurvatureQuadraticField(dt*i, initCurvature);
            } else {
                curvatureTheoretical[i] = Phi.getReferenceCurvatureExplicitEuler(dt, i, initCurvature,  expAngle/180*M_PI, initCP);
            }

            curvatureFile << std::to_string(i*dt) + ", "
                    + std::to_string(curvatureActualDivergence[i]) + ", "
                    + std::to_string(curvatureActualHeight[i]) + ", "
                    + std::to_string(curvatureTheoretical[i]) + "\n";

            std::cout << "Measured curvature with divergence:     " + std::to_string(curvatureActualDivergence[i]) + "\n";
            std::cout << "Measured curvature with height function:" + std::to_string(curvatureActualHeight[i]) + "\n";
            std::cout << "Reference curvature: " + std::to_string(curvatureTheoretical[i])  << std::endl;
        }
        
        // Calculate numerical flux through all faces of each cell update Phi
        Phi.calculateNextTimestep(dt, i);

        positionFile.flush();
        angleFile.flush();
        if (calculateCurvature)
            curvatureFile.flush();
    }
    // Add closing line to the XMF file
    xmfFile << "</Grid>\n</Domain>\n</Xdmf>" << std::endl;

    std::cout << std::endl;
    std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
    std::cout << "Sum of Phi at end: " << Phi.sumLevelSet() << std::endl;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Total runtime: " << duration.count() << "s" << std::endl;


    return 0;
}
