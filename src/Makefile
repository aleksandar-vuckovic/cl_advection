default: finiteVolumes

sourceFiles = main.cpp Field.cpp LevelSet.cpp VelocityField.cpp velocityFields.cpp vecMath3D.cpp
objects = main.o Field.o LevelSet.o VelocityField.o velocityFields.o vecMath3D.o

main.o : main.cpp
	g++ -c -O3 -Wall main.cpp -o main.o

Field.o : Field.cpp
	g++ -c -O3 -Wall Field.cpp -o Field.o

LevelSet.o : LevelSet.cpp
	g++ -c -O3 -Wall LevelSet.cpp -o LevelSet.o

VelocityField.o : VelocityField.cpp
	g++ -c -O3 -Wall VelocityField.cpp -o VelocityField.o

velocityFields.o : velocityFields.cpp
	g++ -c -O3 -Wall velocityFields.cpp -o velocityFields.o

vecMath3D.o : vecMath3D.cpp
	g++ -c -O3 -Wall vecMath3D.cpp -o vecMath3D.o

finiteVolumes: $(objects)
	g++ $(objects) -O3 -Wall -o finiteVolumes

clean:
	rm -r data/ $(objects) finiteVolumes finiteVolumes.debug *csv* *xmf*

debug: $(sourceFiles)
	g++ -g $(sourceFiles) -o finiteVolumes.debug
run:
	./finiteVolumes

brun:
	make
	./finiteVolumes
