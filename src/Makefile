default: finiteVolumes

finiteVolumes.o : finiteVolumes.cpp
	g++ -c -Ofast -Wall finiteVolumes.cpp -o finiteVolumes.o

vecMath3D.o : vecMath3D.cpp
	g++ -c -Ofast -Wall vecMath3D.cpp -o vecMath3D.o

fields.o : fields.cpp
	g++ -c -Ofast -Wall fields.cpp -o fields.o	
	
VelocityField.o : VelocityField.cpp
	g++ -c -Ofast -Wall VelocityField.cpp -o VelocityField.o

finiteVolumes: finiteVolumes.o vecMath3D.o fields.o VelocityField.o
	g++ finiteVolumes.o fields.o vecMath3D.o VelocityField.o -Ofast -Wall -o finiteVolumes

clean:
	rm -r data/ finiteVolumes.o vecMath2D.o vecMath3D.o fields.o VelocityField.o finiteVolumes finiteVolumes.debug *csv* *xmf*

debug: finiteVolumes.cpp vecMath2D.cpp 
	g++ -g finiteVolumes.cpp vecMath3D.cpp fields.cpp VelocityField.cpp -o finiteVolumes.debug
run:
	./finiteVolumes

brun:
	make
	./finiteVolumes
