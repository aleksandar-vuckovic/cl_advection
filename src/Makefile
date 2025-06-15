default: cl_advection

sourceFiles = main.cpp Field.cpp LevelSet.cpp VelocityField.cpp velocityFields.cpp vecMath.cpp
objects = main.o Field.o LevelSet.o VelocityField.o velocityFields.o vecMath.o
objectFlags = -c -O3 -Wall -fopenmp


main.o : main.cpp
	$(CXX) $(objectFlags) -I $(EIGEN_PATH) main.cpp -o main.o

Field.o : Field.cpp
	$(CXX) $(objectFlags) Field.cpp -o Field.o

LevelSet.o : LevelSet.cpp
	$(CXX) $(objectFlags) -I $(EIGEN_PATH) LevelSet.cpp -o LevelSet.o

VelocityField.o : VelocityField.cpp
	$(CXX) $(objectFlags) VelocityField.cpp -o VelocityField.o

velocityFields.o : velocityFields.cpp
	$(CXX) $(objectFlags) velocityFields.cpp -o velocityFields.o

vecMath.o : vecMath.cpp
	$(CXX) $(objectFlags) vecMath.cpp -o vecMath.o

cl_advection: $(objects)
	$(CXX) $(objects) -O3 -Wall -fopenmp -o cl_advection

clean:
	rm -rf data/ $(objects) cl_advection cl_advection.debug cl_advection.profile *csv* *xmf*

debug: $(sourceFiles)
	$(CXX) -I $(EIGEN_PATH) -g -Wall -fopenmp $(sourceFiles) -o cl_advection.debug

profile: $(sourceFiles)
	$(CXX) -pg -Wall -fopenmp $(sourceFiles) -o cl_advection.profile

