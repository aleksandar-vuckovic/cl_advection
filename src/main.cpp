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
#include <getopt.h>   // GNU getopt
#include "LevelSet.hpp"
#include "VelocityField.hpp"
#include "enums.hpp"

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
int main(int argc, char **argv) {
	auto start = std::chrono::system_clock::now();

    int numX, numY, numZ, threads;
    numX = numY = numZ = 0;
    threads = 1;
    
    double lenX, lenY, lenZ, time, centerX, centerY, centerZ, radius, expcpX, expcpY, expcpZ, expAngle, expNormalX, expNormalY, expNormalZ, initCurvature;
    double v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6, tau, CFL, writestepsFraction, planePolarAngle, planeAzimuthalAngle, fieldAzimuthalAngle, alpha;
    double paraboloidStretchX, paraboloidStretchZ, paraboloidHeightMinimum;
    double ellipsoidStretchX, ellipsoidStretchY, ellipsoidStretchZ;

    lenX = lenY = lenZ = time = centerX = centerY = centerZ = radius = expcpX = expcpY = expcpZ = expNormalX = expNormalY = expNormalZ = initCurvature
    = expAngle = v0 = w0 = x0 = y0 = z0 = c1 = c2 = c3 = c4 = c5 = c6 = tau = CFL = writestepsFraction = planePolarAngle
    = planeAzimuthalAngle =  fieldAzimuthalAngle = alpha = paraboloidStretchX = paraboloidStretchZ = paraboloidHeightMinimum 
    = ellipsoidStretchX = ellipsoidStretchY = ellipsoidStretchZ = 0;
    
    bool writeField = false, calculateCurvature = true;
    std::string trackedContactPoint = "left", fieldName = "", outputDirectory = "";
    InitShape initShape = InitShape::sphere;
    VelocityField *field = nullptr;

    // Read data from Inputfile
    std::ifstream inFileStream("Inputfile");
    if(!inFileStream.good()){
        std::cout << "Error: File 'Inputfile' not found in current folder.\n";
        exit(-1);
    }

    std::string line, varName, value, endMessage;

    try {
        while(std::getline(inFileStream, line)) {
            std::istringstream linestream(line);
            if(std::getline(linestream, varName, '=')) {
                // In case a line is commented out
                if (varName.substr(0, 2) == "//")
                    continue;

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
                    else if (varName == "field") {
                        fieldName = value;
                        if (fieldName == "strawberryField" && fieldAzimuthalAngle != 0)
                            throw std::invalid_argument("Error: field is strawberryField and fieldAzimuthalAngle != 0. "
                                                        "It is currently not possible to rotate the strawberry field.");
                    }
                    else if (varName == "fieldAzimuthalAngle") {
                        fieldAzimuthalAngle = std::stod(value);
                        if (fieldName == "strawberryField" && fieldAzimuthalAngle != 0)
                            throw std::invalid_argument("Error: field is strawberryField and fieldAzimuthalAngle != 0. "
                                                        "It is currently not possible to rotate the strawberry field.");
                    }
                    else if (varName == "alpha")
                        alpha = std::stod(value);
                    else if (varName == "geometryType" || varName == "initShape") {
                        if (value == "sphere")
                            initShape = InitShape::sphere;
                        else if (value == "plane")
                            initShape = InitShape::plane;
                        else if (value == "paraboloid")
                            initShape = InitShape::paraboloid;
                        else if (value == "ellipsoid")
                            initShape = InitShape::ellipsoid;
                        else
                            throw std::invalid_argument("No valid initialization shape chosen. Please choose either sphere, plane, paraboloid or ellipsoid.");
                    }
                    else if (varName == "centerX")
                        centerX = std::stod(value);
                    else if (varName == "centerY")
                        centerY = std::stod(value);
                    else if (varName == "centerZ")
                        centerZ = std::stod(value);
                    else if (varName == "radius") {

                        if (initShape == InitShape::sphere) {
                            radius = std::stod(value);
                        } else {
                            throw std::invalid_argument("Given radius while not initializing sphere.");
                        }
                    }
                    else if (varName == "planePolarAngle") {
                        if (initShape == InitShape::plane) {
                            planePolarAngle = std::stod(value);
                        } else {
                            throw std::invalid_argument("Given plane angle while not initializing plane.");
                        }
                    }
                    else if (varName == "planeAzimuthalAngle") {
                        if (initShape == InitShape::plane) {
                            planeAzimuthalAngle = std::stod(value);
                        } else {
                            throw std::invalid_argument("Given plane angle while not initializing plane.");
                        }
                    }
                    else if (varName == "paraboloidStretchX")
                        paraboloidStretchX = std::stod(value);
                    else if (varName == "paraboloidStretchZ")
                        paraboloidStretchZ = std::stod(value);
                    else if (varName == "paraboloidHeightMinimum")
                        paraboloidHeightMinimum = std::stod(value);
                    else if (varName == "ellipsoidStretchX")
                        ellipsoidStretchX = std::stod(value);
                    else if (varName == "ellipsoidStretchY")
                        ellipsoidStretchY = std::stod(value);
                    else if (varName == "ellipsoidStretchZ")
                        ellipsoidStretchZ = std::stod(value);
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
                    else if (varName == "expNormalX" || varName == "expNormalY" || varName == "expNormalZ") {
                        endMessage = "WARNING: Parameters expNormalX / expNormalY / expNormalZ are deprecated.\n In 3D, the expected normal vector "
                                      "is instead calculated from the given contact point and the shape arguments. In 2D, expAngle is still used.";
                    }
                    else if (varName == "trackedContactPoint")
                        trackedContactPoint = value;
                    else if (varName == "initCurvature")
                        initCurvature = std::stod(value);
                    else if (varName == "calculateCurvature")
                        std::stringstream(value) >> std::boolalpha >> calculateCurvature;
                    else if (varName == "threads")
                        threads = std::stoi(value);
                    else
                        throw std::invalid_argument("Input parameter \"" + varName + "\" not recognized.");
                    }
                else
                    throw std::invalid_argument("Input parameter \"" + varName + "\" is not set.");
            }
        }
    } catch (std::invalid_argument &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    static struct option long_options[] = {
        {"resolution", required_argument, nullptr, 'r'},
        {"field", required_argument, nullptr, 'f'},
        {"CFL",required_argument , nullptr, 'c'},
        {"output", required_argument, nullptr, 'o'},
        {"threads", required_argument, nullptr, 't'},
        {"writeField", required_argument, nullptr, 'w'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    extern char* optarg;
    // Parse command line arguments with GNU getopt
    while (true) {
        int option_index;
        opt = getopt_long(argc, argv, "f:r:c:o:t:w:", long_options, &option_index);

        if (opt == -1)      //If there are no more options or arguments left, exit the loop
            break;

        switch(opt) {
            case 'f': {
                fieldName = optarg;
                break;
            }
            case 'r': {
                double resolution = std::stod(optarg);
                numX = numX * resolution;
                numY = numY * resolution;
                if (numZ != 1)
                    numZ = numZ * resolution;
                break;
            }
            case 'c': {
                CFL = std::stod(optarg);
                break;
            }
            case 'o': {
                outputDirectory = optarg;
                if (outputDirectory.back() != '/')
                    outputDirectory = outputDirectory + '/';
                break;
            }
            case 't': {
                threads = std::stoi(optarg);
                break;
            }

            case 'w': {
                std::stringstream(optarg) >> std::boolalpha >> writeField;
		break;
            }

            default:
                abort();

        }
    }

    while (optind < argc) {
        std::cout << "non-option argument" << std::string(argv[optind]) << std::endl;
        optind++;
    }

  	//Parallel computing
  	omp_set_num_threads(threads);

    field = new VelocityField(fieldName, v0, w0, x0, y0, z0, c1, c2, c3, c4, c5, c6,
                              tau, 0, lenX, 0, lenY, 0, lenZ, lenX/numX,lenY/numY,lenZ/numZ, fieldAzimuthalAngle/180*M_PI, alpha, outputDirectory);

    double dx = lenX/numX;
    double dy = lenY/numY;
    double dz = lenZ/numZ;

    std::vector<double> shapeParams;  // This will be populated below, depending on the chosen inital shape and used in LevelSet::expectedNormalVector

    try {
        if (initShape == InitShape::sphere) {
            if (initCurvature != 0)
                throw std::invalid_argument("Given initCurvature when initializing a sphere. This is not necessary as it can be calculated form the radius.");
            if (numZ == 1) {
                initCurvature = -1/radius;
            } else {
                initCurvature = -2/radius;
            }
            shapeParams = {};    // No additional parameters needed for the sphere
        } else if (initShape == InitShape::plane) {
            if (initCurvature != 0)
                throw std::invalid_argument("Given non-zero initCurvature when initializing a plane.");
            initCurvature = 0;
            shapeParams = {planePolarAngle / 180 * M_PI, planeAzimuthalAngle / 180 * M_PI};
        } else if (initShape == InitShape::paraboloid) {
            shapeParams = {paraboloidStretchX, paraboloidStretchZ};
        } else if (initShape == InitShape::ellipsoid) {
            shapeParams = {ellipsoidStretchX, ellipsoidStretchY, ellipsoidStretchZ};
        }
    }  catch (std::invalid_argument& e) {
        std::cout << e.what() << std::endl;
    }

    Vector center = {centerX, centerY, centerZ};
    Vector expCP = {expcpX, expcpY, expcpZ};

    Matrix expNormalVecGrad; // This will be used to calculate the curvature derivative at t = 0

    double dt = 0.0;
    if ( numZ == 1 && expNormalX == 0 && expNormalY == 0 && expNormalZ == 0 ) {
      // 2D case
      dt = CFL*std::min(dx,dy)/field->getMaxNormValue();
    }
    else {
      // 3D case
      dt = CFL*std::min(std::min(dx,dy),dz)/field->getMaxNormValue();
    } 

    int timesteps = time/dt;
    int writesteps = ceil(writestepsFraction*timesteps);

    LevelSet Phi(numX, numY, numZ, dx, dy, dz, field, trackedContactPoint, dt, timesteps, 
    expCP, expAngle, initCurvature, initShape, shapeParams, center, outputDirectory);

    expNormalVecGrad = Phi.expectedNormalVectorGradient(expCP);

    std::vector<Vector> positionTheoretical = Phi.getPositionReference();
    std::vector<double> angleTheoretical = Phi.getAngleReference();
    std::vector<double> curvatureTheoretical = Phi.getCurvatureReference();
    std::vector<Vector> tangentAReference = Phi.getTangentAReference();
    std::vector<Vector> tangentBReference = Phi.getTangentBReference();
    std::vector<Vector> tangentCReference = Phi.getTangentCReference();
    std::vector<double> sectionalCurvatureAReference = Phi.getSectionalCurvatureAReference();
    std::vector<double> sectionalCurvatureBReference = Phi.getSectionalCurvatureBReference();
    std::vector<double> sectionalCurvatureCReference = Phi.getSectionalCurvatureCReference();
    std::vector<double> angleActual(timesteps);
    std::vector<double> curvatureActualDivergence(timesteps);
    std::vector<double> sectionalCurvatureA_Actual(timesteps);
    std::vector<double> sectionalCurvatureB_Actual(timesteps);
    std::vector<double> sectionalCurvatureC_Actual(timesteps);

    if (initShape == InitShape::sphere)
        Phi.initSphere(center, radius);
    else if (initShape == InitShape::plane)
        Phi.initPlane(center, planePolarAngle, planeAzimuthalAngle);
    else if (initShape == InitShape::paraboloid)
        Phi.initParaboloid(center, paraboloidStretchX, paraboloidStretchZ, paraboloidHeightMinimum);
    else if (initShape == InitShape::ellipsoid)
        Phi.initEllipsoid(center, ellipsoidStretchX, ellipsoidStretchY, ellipsoidStretchZ);

    std::string temp = "mkdir -p " + outputDirectory;
    int sysRet;
    if (outputDirectory != "")
        sysRet = system(temp.data());
    temp = "mkdir " + outputDirectory + "data";
    sysRet = system(temp.data());
    if (sysRet == 256) {
    	std::cout << "Overwriting folder \"data\".\n";
    }
    
    std::ofstream positionFile(outputDirectory + "position.csv");
    positionFile.precision(16);
    std::ofstream angleFile(outputDirectory + "contactAngle.csv");
    angleFile.precision(16);
    std::ofstream curvatureFile(outputDirectory + "curvature.csv");
    curvatureFile.precision(16);
    std::ofstream curvatureDerivativeFile(outputDirectory + "curvatureDerivative.csv");
    curvatureDerivativeFile.precision(16);
    std::ofstream sectionalCurvatureAFile(outputDirectory + "sectionalCurvatureA.csv");
    sectionalCurvatureAFile.precision(16);
    std::ofstream sectionalCurvatureBFile(outputDirectory + "sectionalCurvatureB.csv");
    sectionalCurvatureBFile.precision(16);
    std::ofstream sectionalCurvatureCFile(outputDirectory + "sectionalCurvatureC.csv");
    sectionalCurvatureCFile.precision(16);
    std::ofstream minimumGradient(outputDirectory + "minimumGradient.csv");
    minimumGradient.precision(16);

    double sumAtStart = Phi.sumLevelSet();

    // XMF file for Paraview
    std::ofstream MainXmfFile(outputDirectory + "data/Phi.xmf");
    std::ofstream tauXmfFile(outputDirectory + "data/Tau.xmf");



    // main loop
    for (int i = 0; i < timesteps; ++i) {
        std::cout << "Step " << i << std::endl;
        //Write field to file
        if (writeField && i % (int)ceil((double)timesteps/writesteps) == 0) {
            Phi.writeToFile(dt, i, timesteps, writesteps, &MainXmfFile, &tauXmfFile);
        }

	    // Calculate the reference position of the contact point
        Vector newCPReference = positionTheoretical[i];

        // Get the the new contact point
        Vector newCPActual;
        try {
            newCPActual = Phi.getContactPoint(i);
        } catch (std::runtime_error& e) {
            std::string str = e.what();
            if (str.compare("CP_Exit_Simulation_Plane") == 0) {
                std::cout << "WARNING: Contact point left simulation plane during evolution time. Exiting.\n";
                positionFile.flush();
                angleFile.flush();
                if (calculateCurvature)
                    curvatureFile.flush();
                break;
            }
        }

        // Evaluate the Contact Angle numerically based on Phi
        angleActual[i] = Phi.getContactAngleInterpolated(i);

        // Output to command line and positionFile
        std::cout << "Time: " << std::to_string(i*dt) << "\n";

        if (numZ == 1) {
            positionFile << i*dt << ", " << newCPActual[0] << ", " << newCPReference[0] << std::endl;
        } else {
            positionFile << i*dt << ", "
                         << newCPReference[0] << ", " << newCPReference[1] << ", " << newCPReference[2] << std::endl;
        }
        std::cout << "Position " << "x, z  = " << newCPReference[0] << ", " << newCPReference[2] << "\n";

        std::cout << "Actual: " << std::to_string(angleActual[i]/(2*M_PI)*360) + "\n";
		angleFile << std::to_string(i*dt) << ", "
				  << angleActual[i]/(2*M_PI)*360 << ", "
				  << angleTheoretical[i] << "\n";
		std::cout << "Reference: " << angleTheoretical[i] << "\n"; 
                
        minimumGradient << i*dt << ", " << Phi.getMinimalGradientNorm() << "\n";

        if (calculateCurvature) {

            curvatureActualDivergence[i] = Phi.getCurvatureInterpolated(i);
            sectionalCurvatureA_Actual[i] = Phi.getSectionalCurvatureInterpolated(i, tangentAReference[i]);
            sectionalCurvatureB_Actual[i] = Phi.getSectionalCurvatureInterpolated(i, tangentBReference[i]);
            sectionalCurvatureC_Actual[i] = Phi.getSectionalCurvatureInterpolated(i, tangentCReference[i]);

            curvatureFile << std::to_string(i*dt) + ", "
                    + std::to_string(curvatureActualDivergence[i]) + ", "
                    + std::to_string(curvatureTheoretical[i]) + "\n";

            std::cout << "Measured curvature with divergence: " + std::to_string(curvatureActualDivergence[i]) + "\n";
            std::cout << "Reference curvature: " + std::to_string(curvatureTheoretical[i])  << std::endl;
            std::cout << "Reference secCurvA: " + std::to_string(sectionalCurvatureAReference[i])  << std::endl;

            sectionalCurvatureAFile << std::to_string(i*dt) + ", "
                             + std::to_string(sectionalCurvatureA_Actual[i]) + ", "
                             + std::to_string(sectionalCurvatureAReference[i]) + "\n";

            sectionalCurvatureBFile << std::to_string(i*dt) + ", "
                                       + std::to_string(sectionalCurvatureB_Actual[i]) + ", "
                                       + std::to_string(sectionalCurvatureBReference[i]) + "\n";
            sectionalCurvatureCFile << std::to_string(i*dt) + ", "
                                       + std::to_string(sectionalCurvatureC_Actual[i]) + ", "
                                       + std::to_string(sectionalCurvatureCReference[i]) + "\n";
        }
        
        // Calculate numerical flux through all faces of each cell and update Phi
        Phi.calculateNextTimestep(dt, i);

        positionFile.flush();
        angleFile.flush();
        curvatureFile.flush();
        curvatureDerivativeFile.flush();
        sectionalCurvatureAFile.flush();
        sectionalCurvatureBFile.flush();
        minimumGradient.flush();
    }

    if (numZ > 1) {
        curvatureDerivativeFile << "0, "
                                << std::to_string((curvatureActualDivergence[1] - curvatureActualDivergence[0]) / dt) + ", "
                                << std::to_string(Phi.referenceCurvatureDeriv3D(initCurvature, expNormalVecGrad)) + "\n";
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
    std::cout << endMessage;


    return 0;
}
