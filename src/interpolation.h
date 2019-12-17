#pragma once

#include <sstream>
#include <vector>
#include <iostream>
#include "array"

//the bilinear interpolation class that can be accessed outside of this file
class BilinearInterpolation {

//do not allow these functions and variables to be externally accessible
private:
//    create a static function that returns a pair of index values of the nearest neighbor values
    static std::pair<unsigned int, unsigned int> get_neighbors(
            double val1, //val1 is the value which will get a nearest neighbor from var1_vector
            double val2, //val1 is the value which will get a nearest neighbor from var2_vector
            const std::vector<double> &var1_vector, //the vector from with var1 will be interpolated
            const std::vector<double> &var2_vector //the vector from with var2 will be interpolated
    ) {
        auto line = std::lower_bound(var1_vector.begin(), var1_vector.end(),
                                     val1); //returns an iterator that points to the first value in var1_vector that is not < val1
        const unsigned int line_index = std::distance(var1_vector.begin(),
                                                      line); //returns an int which is the number of elements between the first element of var1_vector and line
        auto column = std::lower_bound(var2_vector.begin(), var2_vector.end(),
                                       val2); //returns an iterator that points to the first value in var2_vector that is not < val2
        const unsigned int column_index = std::distance(var2_vector.begin(),
                                                        column); //returns an int which is the number of elements between the first element of var2_vector and column

        const unsigned int max_line_index = var1_vector.size() - 1; //get the largest index value possible of var1_vector
        const unsigned int max_column_index = var2_vector.size() - 1; //get the largest index value possible of var2_vector

//        make a tuple of the nearest neighbor indices: std::min returns either the value in argument 1 (the index calculated above) or the value of argment 2, which is the max index (i.e. the boundary)
        return std::make_pair(std::min(line_index, max_line_index), std::min(column_index, max_column_index));
    }

//allow these functions and variables to be externally accessible
public:
    static double interpolate(
            double val1,
            double val2,
            const std::vector<double> &var1_vector,
            const std::vector<double> &var2_vector,
            const unsigned int property_index,
            const std::vector<std::vector<std::array<double, 6>>> &eos_data
    ) {
        std::pair<unsigned int, unsigned int> neighbors = get_neighbors(val1, val2, var1_vector, var2_vector);

        // If we are not at the boundaries of the data table
        if ((neighbors.first < var1_vector.size() - 1) &&
            (neighbors.second < var2_vector.size() - 1) &&
            (neighbors.first > 0) &&
            (neighbors.second > 0)) {
            // compute the interpolation weights of this val1 and val2
            const auto xi = (val1 - var1_vector[neighbors.first - 1]) /
                            (var1_vector[neighbors.first] - var1_vector[neighbors.first - 1]);

            const auto eta = (val2 - var2_vector[neighbors.second-1]) /
                             (var2_vector[neighbors.second] - var2_vector[neighbors.second-1]);

            // use these coordinates for a bilinear interpolation
            const double interpolated_value = (1 - xi) * (1 - eta) * eos_data[neighbors.first - 1][neighbors.second - 1][property_index] +
                                              xi * (1 - eta) * eos_data[neighbors.first][neighbors.second - 1][property_index] +
                                              (1 - xi) * eta * eos_data[neighbors.first - 1][neighbors.second][property_index] +
                                              xi * eta * eos_data[neighbors.first][neighbors.second][property_index];
            return interpolated_value;
        } else {
            // Return the boundary value
            const double interpolated_value = eos_data[neighbors.first][neighbors.second][property_index];
            return interpolated_value;
        }
    };
};