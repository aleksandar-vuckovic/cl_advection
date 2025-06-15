default: all

shear_2D: shear_2D.cpp helpers.hpp
	$(CXX) -g shear_2D.cpp helpers.hpp -o shear_2D

navier_2D: navier_2D.cpp helpers.hpp
	$(CXX) -g navier_2D.cpp helpers.hpp -o navier_2D


shear_3D: shear_3D.cpp helpers.hpp
	$(CXX) -g shear_3D.cpp helpers.hpp -o shear_3D

strawberry_3D: strawberry_3D.cpp helpers.hpp
	$(CXX) -g strawberry_3D.cpp helpers.hpp -o strawberry_3D

all: shear_2D navier_2D shear_3D strawberry_3D

clean:
	rm ./shear_2D ./navier_2D ./shear_3D ./strawberry_3D
