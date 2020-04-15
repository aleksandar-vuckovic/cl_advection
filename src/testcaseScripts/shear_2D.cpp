#include <iostream>
#include "helpers.hpp"

int main() {

  std::string path = "../testcases/2D/curvature/convergenceAnalysis/shear/";
  Data curvatureData_25 = read_csv(path + "shear_25/curvature.csv");
 Data curvatureData_50 = read_csv(path + "shear_50/curvature.csv");
  Data curvatureData_100 = read_csv(path + "shear_100/curvature.csv");
  Data curvatureData_200 = read_csv(path + "shear_200/curvature.csv");

  Data angleData_25 = read_csv(path + "shear_25/contactAngle.csv");
  Data angleData_50 = read_csv(path + "shear_50/contactAngle.csv");
  Data angleData_100 = read_csv(path + "shear_100/contactAngle.csv");
  Data angleData_200 = read_csv(path + "shear_200/contactAngle.csv");


  double time_tolerance = 0.1;
  
  double value_c = 0.15;
  double value_c_tolerance = 0.01;
  
  bool success_curvature_25 = verify(curvatureData_25, 5, time_tolerance, value_c, value_c_tolerance);
  bool success_curvature_50 = verify(curvatureData_50, 5, time_tolerance, value_c/2, value_c_tolerance);
  bool success_curvature_100 = verify(curvatureData_100, 5, time_tolerance, value_c/4, value_c_tolerance);
  bool success_curvature_200 = verify(curvatureData_200, 5, time_tolerance, value_c/8, value_c_tolerance);

  double value_a = 0.6;
  double value_a_tolerance = 0.1;
  
  bool success_angle_25 = verify(angleData_25, 3, time_tolerance, value_a, value_a_tolerance);
  bool success_angle_50 = verify(angleData_50, 3, time_tolerance, value_a/2, value_a_tolerance);
  bool success_angle_100 = verify(angleData_100, 3, time_tolerance, value_a/4, value_a_tolerance);
  bool success_angle_200 = verify(angleData_200, 3, time_tolerance, value_a/8, value_a_tolerance);

  if (success_curvature_25 && success_curvature_50 && success_curvature_100 && success_curvature_200
      && success_angle_25 && success_angle_50 && success_angle_100 && success_angle_200) {
    
    std::cout << "Test passed!\n";
    return 0;
    
  } else {
    std::cout << "Test failed.\n";
    return 1;
  }
}

