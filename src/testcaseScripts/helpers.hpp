#ifndef HELPERS_HPP
#define HELPERS_CPP

#include <fstream>
#include <vector>
#include <sstream>



typedef std::vector<std::vector<double>> Data;

Data read_csv(std::string path) {

  std::ifstream stream(path);
  std::string line;
  Data tempReturn;

  while(std::getline(stream, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    std::vector<double> readLine;
    
    while(std::getline(lineStream, cell, ',')) {
      readLine.push_back(std::stod(cell));
    }

    tempReturn.push_back(readLine);
  }

  return tempReturn;
}

bool verify(Data data, double time, double time_tolerance, double exp_value, double value_tolerance) {

  for (auto& row: data) {
    double _time_ = row[0];
    if (std::abs(_time_ - time) < time_tolerance) {
      double value = std::abs(row[1] - row[2]);

      return std::abs(value - exp_value) < value_tolerance;
    }
  }
  throw std::invalid_argument("No point with given time value found.");
}

#endif
