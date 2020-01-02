#pragma once

#include <sstream>
#include <vector>
#include <iostream>
#include "array"

template<typename T>
void print(std::vector<T> const &v)
{
    for (auto i: v) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
}

template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;

    std::vector<T> vec(first, last);
    return vec;
}

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

        const unsigned int max_line_index = var1_vector.size() - 1; // get the largest index value possible of var1_vector
        const unsigned int max_column_index = var2_vector.size() - 1; //get the largest index value possible of var2_vector

//        make a tuple of the nearest neighbor indices: std::min returns either the value in argument 1 (the index calculated above) or the value of argument 2, which is the max index (i.e. the boundary)
        return std::make_pair(std::min(line_index, max_line_index), std::min(column_index, max_column_index));
    }

    static std::pair<unsigned int, unsigned int> get_neighbors_restricted_indices(
            double val, //val1 is the value which will get a nearest neighbor from var1_vector
            const std::vector<double> &var_vector //the vector from with var1 will be interpolated
            ) {

        auto lower_bound = std::lower_bound(var_vector.begin(), var_vector.end(),
                                     val); //returns an iterator that points to the first value in var1_vector that is not < val
        const unsigned int lower_bound_index = std::distance(var_vector.begin(), lower_bound);
        auto upper_bound = std::upper_bound(var_vector.begin(), var_vector.end(),
                                            val);
        const unsigned int upper_bound_index = std::distance(var_vector.begin(), upper_bound);

        return std::make_pair(lower_bound_index, upper_bound_index);
    }


//allow these functions and variables to be externally accessible
public:
    static double interpolate(
            double val1,
            double val2,
            const std::vector<double> &var1_vector,
            const std::vector<double> &var2_vector,
            const unsigned int property_index,
            const std::vector<std::vector<std::array<double, 6>>> &eos_data,
            bool restrict
    ) {
        std::pair<unsigned int, unsigned int> neighbors;

        // if we are *not* restricting indices, then just get neighbors within the entire vector
        if (!restrict) {
            neighbors = get_neighbors(
                    val1,
                    val2,
                    var1_vector,
                    var2_vector
                    );
        }
        // if we *are* restricting indices, then just get neighbors within the restricted vector
        else {
            std::pair<unsigned int, unsigned int> restricted_index_pair = get_neighbors_restricted_indices(var1_vector, 1.02232811E+02, 120);
//            std::cout << "************" << std::endl;
//            std::cout << restricted_index_pair.first << std::endl;
//            std::cout << restricted_index_pair.second << std::endl;
//            std::cout << "************" << std::endl;
            neighbors = get_neighbors(
                    val1,
                    val2,
                    slice(var1_vector, restricted_index_pair.first, restricted_index_pair.second),
                    slice(var2_vector, restricted_index_pair.first, restricted_index_pair.second)
                    );
        }

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

    static std::pair<unsigned int, unsigned int> get_neighbors_restricted_indices(
            const std::vector<double> &array,
            double value,
            unsigned int grid_length
    ) {
        std::cout << array.size() << std::endl;
        // this function returns the upper and lower indices to restrict an array to 3 values: the nearest value to value, the value right above value, and the value right below value
        // assume that the array is ordered

        unsigned int lower_bound_index = 0;  // instantiate the lower bound
        unsigned int upper_bound_index = 0;  // instantiate the upper bound
        std::vector<double> reversed_vector(array.size());
        std::reverse_copy(std::begin(array), std::end(array), std::begin(reversed_vector)); // reverse the array

//        std::cout << "myvector contains:" << std::endl;
//        for (const auto& v : reversed_vector) {
//            std::cout << v << " " << std::endl;
//        }

        // search for the upper bound index value
        for (unsigned int a = 0; a < array.size(); a = a + 1) {
            double v = array[a];
            if (v > value) {
                upper_bound_index = a + grid_length - 1;
                if (upper_bound_index > array.size()) {
                    upper_bound_index = array.size() - 1;
                }
                break;
            }
        }

        // search for the lower bound index value
        for (unsigned int a = 0; a < reversed_vector.size(); a = a + 1) {
            double v = reversed_vector[a];
            if (v < value) {
                lower_bound_index = array.size() - (a - 1) - grid_length;  // fetch the index value from the original array, not the reversed array
                if (lower_bound_index < 0) {
                    lower_bound_index = 0;
                }
                break;
            }
        }

        // return a pair of the lower and upper bound index values
        return std::make_pair(lower_bound_index, upper_bound_index);
    }
};