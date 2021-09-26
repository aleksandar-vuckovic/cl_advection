default: cl_advection

sourceFiles = main.cpp Field.cpp LevelSet.cpp VelocityField.cpp velocityFields.cpp vecMath.cpp
objects = main.o Field.o LevelSet.o VelocityField.o velocityFields.o vecMath.o
objectFlags = -c -O3 -Wall -fopenmp

main.o : main.cpp
	g++ $(objectFlags) -I /usr/include/eigen3 main.cpp -o main.o

Field.o : Field.cpp
	g++ $(objectFlags) Field.cpp -o Field.o

LevelSet.o : LevelSet.cpp
	g++ $(objectFlags) -I /usr/include/eigen3 LevelSet.cpp -o LevelSet.o

VelocityField.o : VelocityField.cpp
	g++ $(objectFlags) VelocityField.cpp -o VelocityField.o

velocityFields.o : velocityFields.cpp
	g++ $(objectFlags) velocityFields.cpp -o velocityFields.o

vecMath.o : vecMath.cpp
	g++ $(objectFlags) vecMath.cpp -o vecMath.o

cl_advection: $(objects)
	g++ $(objects) -O3 -Wall -fopenmp -o cl_advection

clean:
	rm -rf data/ $(objects) cl_advection cl_advection.debug cl_advection.profile *csv* *xmf*

debug: $(sourceFiles)
	g++ -I /usr/include/eigen3 -g -Wall -fopenmp $(sourceFiles) -o cl_advection.debug

profile: $(sourceFiles)
	g++ -pg -Wall -fopenmp $(sourceFiles) -o cl_advection.profile

