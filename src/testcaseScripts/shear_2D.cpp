#include <iostream>
#include "helpers.hpp"

int main() {

  double time_tolerance = 0.1;
  double value_tolerance = 0.01;
  
  std::string path = "../testcases/2D/curvature/convergenceAnalysis/shear/";
  Data curvatureData_25 = read_csv(path + "shear_25/curvature.csv");
  Data curvatureData_50 = read_csv(path + "shear_50/curvature.csv");
  Data curvatureData_100 = read_csv(path + "shear_100/curvature.csv");
  Data curvatureData_200 = read_csv(path + "shear_200/curvature.csv");

  Data angleData_25 = read_csv(path + "shear_25/contactAnglecsv");
  Data angleData_50 = read_csv(path + "shear_50/contactAnglecsv");
  Data angleData_100 = read_csv(path + "shear_100/contactAnglecsv");
  Data angleData_200 = read_csv(path + "shear_200/contactAnglecsv");

  bool success_curvature_25 = verify(curvatureData_25, 5, time_tolerance, 0.15, value_tolerance);

  if (success_curvature_25)
    std::cout << "Test passed!\n";
  else
    std::cout << "Test failed.\n";
}

