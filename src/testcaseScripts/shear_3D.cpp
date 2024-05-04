#include <iostream>
#include "helpers.hpp"

int main() {

  std::string path = "../testcases/3D/curvature/convergenceAnalysis/shear/";
  Data curvatureData_25 = read_csv(path + "shear_25/curvature.csv");
  Data curvatureData_50 = read_csv(path + "shear_50/curvature.csv");
  Data curvatureData_100 = read_csv(path + "shear_100/curvature.csv");
  Data curvatureData_200 = read_csv(path + "shear_200/curvature.csv");

  Data angleData_25 = read_csv(path + "shear_25/contactAngle.csv");
  Data angleData_50 = read_csv(path + "shear_50/contactAngle.csv");
  Data angleData_100 = read_csv(path + "shear_100/contactAngle.csv");
  Data angleData_200 = read_csv(path + "shear_200/contactAngle.csv");


  double time_tolerance = 0.05;

  double value_a = 1.5;
  double value_a_tolerance = 0.12;
  
  bool success_angle_25 = verify(angleData_25, 0.8, time_tolerance, value_a, value_a_tolerance);
  bool success_angle_50 = verify(angleData_50, 0.8, time_tolerance, value_a/2, value_a_tolerance);
  bool success_angle_100 = verify(angleData_100, 0.8, time_tolerance, value_a/4, value_a_tolerance);
  bool success_angle_200 = verify(angleData_200, 0.8, time_tolerance, value_a/8, value_a_tolerance);

  if (success_angle_25 && success_angle_50 && success_angle_100 && success_angle_200) {
    
    std::cout << "Test passed!\n";
    return 0;
    
  } else {
    std::cout << "Test failed.\n";
    return 1;
  }
}

