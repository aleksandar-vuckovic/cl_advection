default: all

shear_2D: shear_2D.cpp helpers.hpp
	g++ shear_2D.cpp helpers.hpp -o shear_2D

navier_2D: navier_2D.cpp helpers.hpp
	g++ navier_2D.cpp helpers.hpp -o navier_2D

all: shear_2D navier_2D

clean:
	rm ./shear_2D ./navier_2D
