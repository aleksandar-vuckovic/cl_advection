/**
 * @class Streamlines
 * The Streamlines class.
 *
 * Generates stream lines for a static vector field.
 */

#include "Streamlines.hpp"

using std::array;
using std::vector;

/**
 * The constructor.
 * Calculates the streamline field with Runge-Kutta 4 and saves it to the objects data member.
 * @param field The underlying scalar field, for example the levelset field.
 * @param vel The velocity field to generate the streamlines on.
 * @param dt Timestep width to use for streamline calculation.
 */
Streamlines::Streamlines(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField* field, double dt, std::string outputDirectory) :
    Field<int>(numX, numY, numZ, dx, dy, dz) {
    this->outputDirectory = outputDirectory;
	// Array to store streamline points for 10 streamlines
	array<vector< array<double, 3>> , 10> streamlines;

	//Initialize first streamline points
	for (int i = 0; i < 10; ++i) {
        streamlines[i].push_back({field->getXMax()/2, i*0.1*field->getYMax()});
	}

	for (int i = 0; i < 10; ++i) {
		while (true) {
			array<double, 3> p = streamlines[i].back();
            array<double, 3> k1 = dt * field->at(0, p[0], p[1], 0);
            array<double, 3> k2 = dt * field->at(0, p[0] + k1[0]/2, p[1] + k1[1]/2, 0);
            array<double, 3> k3 = dt * field->at(0, p[0] + k2[0]/2, p[1] + k2[1]/2, 0);
            array<double, 3> k4 = dt * field->at(0, p[0] + k3[0], p[1] + k3[1], 0);

			array<double, 3> p_next = p + 1.0/6 * (k1 + 2*k2 + 2*k3 + k4);
            if (p_next[0] < 0 || p_next[0] >= field->getXMax() || p_next[1] < 0 || p_next[1] >= field->getYMax()
				|| streamlines[i].size() > numX*numY/10.0) {
				break;
			}

			streamlines[i].push_back(p_next);

		}
	}

	for (auto line : streamlines)
		for (auto p : line) {
            int i = p[0] / field->getDx();
            int j = p[1] / field->getDy();
			this->at(i, j, 0) = 1;
		}
}

/**
 * Writes the stream line field to disk.
 *
 * Writes the stream line field for visualization in Paraview. Since no XMF file is written,
 * simply calling this function is not enough for visualization. Thus, this function is called within
 * LevelSet::writeToFile, which does write a XMF file.
 *
 * @param t The time
 */
void Streamlines::writeToFile() {
	int Npoints = numX*numY*numZ;
	FILE *streamfile;
    std::string temp = outputDirectory + "data/stream.bin";
    streamfile = fopen(temp.data(), "wb");
	fwrite(getData().data(), sizeof(int), Npoints, streamfile);
	fclose(streamfile);
}
