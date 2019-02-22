default: finiteVolumes

finiteVolumes.o : finiteVolumes.cpp
	g++ -c finiteVolumes.cpp -o finiteVolumes.o

vecMath2D.o : vecMath2D.cpp
	g++ -c vecMath2D.cpp -o vecMath2D.o

finiteVolumes: finiteVolumes.o vecMath2D.o
	g++ finiteVolumes.o vecMath2D.o -o finiteVolumes

clean:
	rm finiteVolumes.o vecMath2D.o finiteVolumes
run:
	./finiteVolumes

brun:
	make
	./finiteVolumes
