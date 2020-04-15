#include <iostream>
#include "helpers.hpp"

int main() {

  std::string path = "../testcases/2D/curvature/convergenceAnalysis/navier/";
  Data curvatureData_25 = read_csv(path + "navier_25/curvature.csv");
  Data curvatureData_50 = read_csv(path + "navier_50/curvature.csv");
  Data curvatureData_100 = read_csv(path + "navier_100/curvature.csv");
  Data curvatureData_200 = read_csv(path + "navier_200/curvature.csv");


  double time = 1.8;
  double time_tolerance = 0.01;
  
  double value_c = 0.011;
  double value_c_tolerance = 0.001;
  
  bool success_curvature_25 = verify(curvatureData_25, time, time_tolerance, value_c, value_c_tolerance);
  bool success_curvature_50 = verify(curvatureData_50, time, time_tolerance, value_c/2, value_c_tolerance);
  bool success_curvature_100 = verify(curvatureData_100, time, time_tolerance, value_c/4, value_c_tolerance);
  bool success_curvature_200 = verify(curvatureData_200, time, time_tolerance, value_c/8, value_c_tolerance);

  if (success_curvature_25 && success_curvature_50 && success_curvature_100 && success_curvature_200) {
    
    std::cout << "Test passed!\n";
    return 0;
    
  } else {
    std::cout << "Test failed.\n";
    return 1;
  }
}

