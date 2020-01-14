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
    numX = numY = numZ = 0;
    threads = 1;
    
    double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius, expcpX, expcpY, expcpZ, expAngle, expNormalX, expNormalY, expNormalZ;
    double v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6, tau, CFL, writestepsFraction, polarAngle, planeAzimuthalAngle, fieldAzimuthalAngle;
    lenX = lenY = lenZ = time = centerX = centerY = centerZ = radius = expcpX = expcpY = expcpZ = expNormalX = expNormalY = expNormalZ = 0;
    expAngle = v0 = w0 = x0 = y0 = z0 = c1 = c2 = c3 = c4 = c5 = c6 = tau = CFL = writestepsFraction = polarAngle = planeAzimuthalAngle =  fieldAzimuthalAngle = 0;
    
    bool writeField = false, calculateCurvature = false;
    std::string trackedContactPoint = "left", fieldName = "", geometryType = "sphere";
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
                    lenX =  std::stod(value);
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
                else if (varName == "w0")
                    w0 = std::stod(value);
                else if (varName == "c1")
                    c1 = std::stod(value);
                else if (varName == "c2")
                    c2 = std::stod(value);
                else if (varName == "c3")
                    c3 = std::stod(value);
                else if (varName == "c4")
                    c4 = std::stod(value);
                else if (varName == "c5")
                    c5 = std::stod(value);
                else if (varName == "c6")
                    c6= std::stod(value);
                else if (varName == "x0")
                    x0 = std::stod(value);
                else if (varName == "y0")
                    y0 = std::stod(value);
                else if (varName == "z0")
                    z0 = std::stod(value);
                else if (varName == "tau")
                    tau = std::stod(value);
                else if (varName == "field")
                    fieldName = value;
                else if (varName == "fieldAzimuthalAngle")
                    fieldAzimuthalAngle = std::stod(value);
                else if (varName == "geometryType")
                    geometryType = value;
                else if (varName == "centerX")
                    centerX = std::stod(value);
                else if (varName == "centerY")
                    centerY = std::stod(value);
                else if (varName == "centerZ")
                    centerZ = std::stod(value);
                else if (varName == "radius") {
                    if (geometryType == "sphere")
                        radius = std::stod(value);
                    else
                        throw std::invalid_argument("Given radius while initializing plane.");
                }
                else if (varName == "polarAngle") {
                    if (geometryType == "plane")
                        polarAngle = std::stod(value);
                    else
                        throw std::invalid_argument("Given plane angle while initializing sphere.");
                }
                else if (varName == "planeAzimuthalAngle") {
                    if (geometryType == "plane")
                        planeAzimuthalAngle = std::stod(value);
                    else
                        throw std::invalid_argument("Given plane angle while initializing sphere.");
                }
                else if (varName == "expcpX")
                    expcpX = std::stod(value);
                else if (varName == "expcpY")
                    expcpY = std::stod(value);
                else if (varName == "expcpZ")
                    expcpZ = std::stod(value);
                else if (varName == "expAngle") {
                    if (numZ == 1)
                        expAngle = std::stod(value);
                    else
                        throw std::invalid_argument("Input parameter expAngle is not applicable in 3D. Instead define the expected normal vector "
                                                    "by setting the parameters expNormalX, expNormalY and expNormalZ.");
                }
                else if (varName == "expNormalX")
                    expNormalX = std::stod(value);
                else if (varName == "expNormalY")
                    expNormalY = std::stod(value);
                else if (varName == "expNormalZ")
                    expNormalZ = std::stod(value);
                else if (varName == "trackedContactPoint")
                    trackedContactPoint = value;
                else if (varName == "calculateCurvature")
                    std::stringstream(value) >> std::boolalpha >> calculateCurvature;
                else if (varName == "threads")
                    threads = std::stoi(value);
                else
                    throw std::invalid_argument("Input parameter \"" + varName + "\" not recognized");
                }
            // In case a line is commented out
            else if (varName.substr(0, 2) == "//")
                continue;
            else
                throw std::invalid_argument("Input parameter \"" + varName + "\" is not set.");
            }
        }

  	//Parallel computing
  	omp_set_num_threads(threads);

    field = new VelocityField(fieldName, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6,
                              tau, 0, lenX, 0, lenY, 0, lenZ, lenX/numX,lenY/numY,lenZ/numZ, fieldAzimuthalAngle);

    double dx = lenX/numX;
    double dy = lenY/numY;
    double dz = lenZ/numZ;
    double initCurvature = 0;
    if (geometryType == "sphere")
        initCurvature = -1/radius;
    else if (geometryType == "plane")
        initCurvature = 0;
    
    array<double, 3> expNormalVec;

    double dt = 0.0;
    if (numZ == 1) {
      // 2D case
      dt = CFL*std::min(dx,dy)/field->getMaxNormValue();

      if (trackedContactPoint == "left") {
          expNormalVec = { -sin(expAngle/180*M_PI), cos(expAngle/180*M_PI), 0};
      } else {
          expNormalVec = {sin(expAngle/180*M_PI), cos(expAngle/180*M_PI), 0};
      }
    }
    else {
      // 3D case
      dt = CFL*std::min(std::min(dx,dy),dz)/field->getMaxNormValue();
      expNormalVec = {expNormalX, expNormalY, expNormalZ};
      expNormalVec = expNormalVec/abs(expNormalVec);
    } 

    int timesteps = time/dt;
    int writesteps = ceil(writestepsFraction*timesteps);

    array<double, 3> center = {centerX, centerY, centerZ};
    array<double, 3> expCP = {expcpX, expcpY, expcpZ};

    LevelSet Phi(numX, numY, numZ, dx, dy, dz, field, trackedContactPoint, dt, timesteps, expCP, expNormalVec, expAngle, initCurvature);
    Streamlines streamlines(numX, numY, numZ, *field, dt);

    std::vector< array<double, 3>> positionTheoretical = Phi.getPositionReference();
    std::vector<double> angleTheoretical = Phi.getAngleReference();
    std::vector<double> curvatureTheoretical = Phi.getCurvatureReference();
    std::vector<double> angleActual(timesteps);
    std::vector<double> curvatureActualDivergence(timesteps);

    if (geometryType == "sphere")
        Phi.initDroplet(center, radius);
    else if (geometryType == "plane") 
        Phi.initPlane(center, polarAngle, planeAzimuthalAngle);

    int sysRet = system("mkdir data");

    if (sysRet == 256) {
    	std::cout << "Overwriting folder \"data\".\n";
    }
    
    std::ofstream positionFile("position.csv");
    positionFile.precision(16);
    std::ofstream angleFile("contactAngle.csv");
    angleFile.precision(16);
    std::ofstream curvatureFile;
    if (calculateCurvature) {
        curvatureFile = std::ofstream("curvature.csv");
        curvatureFile.precision(16);
    }
    double sumAtStart = Phi.sumLevelSet();

    // XMF file for Paraview
    std::ofstream MainXmfFile("data/Phi.xmf");
    std::ofstream tauXmfFile("data/Tau.xmf");



    // main loop
    for (int i = 0; i < timesteps; ++i) {
        std::cout << "Step " << i << std::endl;
        //Write field to file
        if (writeField && i % (int)ceil((double)timesteps/writesteps) == 0) {
            if (i == 0)
                streamlines.writeToFile();
            Phi.writeToFile(dt, i, timesteps, writesteps, &MainXmfFile, &tauXmfFile);
        }

	// Calculate the reference position of the contact point
        array<double, 3> newCPReference;
        if (numZ == 1 && (field->getName() == "navierField" || field->getName() == "timeDependentNavierField"))
        	newCPReference = Phi.contactPointLinearField(dt*i, c1, expcpX, v0);
        else
            newCPReference = positionTheoretical[i];

        // Get the the new contact point
        array<double, 3> newCPActual = Phi.getContactPoint(i);
        
        // Get the indices of the new contact point
        array<int, 3> CP_indices = Phi.getContactPointIndices(i);

        // Evaluate te Contact Angle numerically based on Phi
        angleActual[i] = Phi.getContactAngle(CP_indices);

        // Output to command line and positionFile
        std::cout << "Time: " << std::to_string(i*dt) << "\n";

        if (numZ == 1) {
            positionFile << i*dt << ", " << newCPActual[0] << ", " << newCPReference[0] << std::endl;
        } else {
            positionFile << i*dt << ", "
                         << newCPReference[0] << ", " << newCPReference[1] << ", " << newCPReference[2] << std::endl;
        }

        std::cout << "Actual: " << std::to_string(angleActual[i]/(2*M_PI)*360) + "\n";
		angleFile << std::to_string(i*dt) << ", "
				  << angleActual[i]/(2*M_PI)*360 << ", "
				  << angleTheoretical[i] << "\n";
		std::cout << "Reference: " << angleTheoretical[i] << "\n"; 
                
        
        if (calculateCurvature) {
            curvatureActualDivergence[i] = Phi.getCurvatureDivergence(CP_indices);

            curvatureFile << std::to_string(i*dt) + ", "
                    + std::to_string(curvatureActualDivergence[i]) + ", "
                    + std::to_string(curvatureTheoretical[i]) + "\n";

            std::cout << "Measured curvature with divergence: " + std::to_string(curvatureActualDivergence[i]) + "\n";
            std::cout << "Reference curvature: " + std::to_string(curvatureTheoretical[i])  << std::endl;
        }
        
        // Calculate numerical flux through all faces of each cell and update Phi
        Phi.calculateNextTimestep(dt, i);

        positionFile.flush();
        angleFile.flush();
        if (calculateCurvature)
            curvatureFile.flush();
    }
    // Add closing line to the XMF files
    MainXmfFile << "</Grid>\n</Domain>\n</Xdmf>" << std::endl;
    tauXmfFile << "</Grid>\n</Domain>\n</Xdmf>" << std::endl;

    std::cout << std::endl;
    std::cout << "Sum of Phi at start: " << sumAtStart << std::endl;
    std::cout << "Sum of Phi at end: " << Phi.sumLevelSet() << std::endl;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Total runtime: " << duration.count() << "s" << std::endl;


    return 0;
}
