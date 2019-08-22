default: cl_advection

sourceFiles = main.cpp Field.cpp LevelSet.cpp VelocityField.cpp velocityFields.cpp vecMath.cpp Streamlines.cpp
objects = main.o Field.o LevelSet.o VelocityField.o velocityFields.o vecMath.o Streamlines.o

main.o : main.cpp
	g++ -c -O3 -Wall main.cpp -o main.o

Field.o : Field.cpp
	g++ -c -O3 -Wall Field.cpp -o Field.o

LevelSet.o : LevelSet.cpp
	g++ -c -O3 -Wall -fopenmp LevelSet.cpp -o LevelSet.o

VelocityField.o : VelocityField.cpp
	g++ -c -O3 -Wall VelocityField.cpp -o VelocityField.o

velocityFields.o : velocityFields.cpp
	g++ -c -O3 -Wall velocityFields.cpp -o velocityFields.o

vecMath.o : vecMath.cpp
	g++ -c -O3 -Wall vecMath.cpp -o vecMath.o

Streamlines.o : Streamlines.cpp
	g++ -c -O3 -Wall Streamlines.cpp -o Streamlines.o

cl_advection: $(objects)
	g++ $(objects) -O3 -Wall -fopenmp -o cl_advection

clean:
	rm -rf data/ $(objects) cl_advection cl_advection.debug *csv* *xmf*

debug: $(sourceFiles)
	g++ -g $(sourceFiles) -o cl_advection.debug

profile: $(sourceFiles)
	g++ -pg -fopenmp $(sourceFiles) -o cl_advection.profile

