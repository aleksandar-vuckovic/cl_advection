default: finiteVolumes

finiteVolumes.o : finiteVolumes.cpp
	g++  -c finiteVolumes.cpp -o finiteVolumes.o

finitevolumes: finitevolumes.o
	g++ finiteVolumes.o -o finiteVolumes

clean:
	rm finiteVolumes.o finiteVolumes
run:
	./finiteVolumes
