default: finiteVolumes

finiteVolumes.o : finiteVolumes.cpp
	g++ -c -Ofast -Wall finiteVolumes.cpp -o finiteVolumes.o

vecMath3D.o : vecMath3D.cpp
	g++ -c -Ofast -Wall vecMath3D.cpp -o vecMath3D.o

finiteVolumes: finiteVolumes.o vecMath3D.o
	g++ finiteVolumes.o vecMath3D.o -Ofast -Wall -o finiteVolumes

clean:
	rm -r data/ finiteVolumes.o vecMath2D.o finiteVolumes finiteVolumes.debug *csv* *xmf*

debug: finiteVolumes.cpp vecMath2D.cpp 
	g++ -g finiteVolumes.cpp vecMath3D.cpp -o finiteVolumes.debug
run:
	./finiteVolumes

brun:
	make
	./finiteVolumes
