#pragma once

#include <sstream>
#include <vector>
#include <iostream>
#include "math.h"
#include "array"

// function to print vector contents
template<typename T>
void print(std::vector<T> const &v) {
    for (auto i: v) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
}

// subsamples a vector "v" given a lower index bound "m" and upper index bound "n"
template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n) {
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

        const unsigned int max_line_index =
                var1_vector.size() - 1; // get the largest index value possible of var1_vector
        const unsigned int max_column_index =
                var2_vector.size() - 1; //get the largest index value possible of var2_vector

//        make a tuple of the nearest neighbor indices: std::min returns either the value in argument 1 (the index calculated above) or the value of argument 2, which is the max index (i.e. the boundary)
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
        std::pair<unsigned int, unsigned int> neighbors;

        // if we are *not* restricting indices, then just get neighbors within the entire vector
        neighbors = get_neighbors(
                val1,
                val2,
                var1_vector,
                var2_vector
        );

        // If we are not at the boundaries of the data table
        if ((neighbors.first < var1_vector.size() - 1) &&
            (neighbors.second < var2_vector.size() - 1) &&
            (neighbors.first > 0) &&
            (neighbors.second > 0)) {
            // compute the interpolation weights of this val1 and val2
            const auto xi = (val1 - var1_vector[neighbors.first - 1]) /
                            (var1_vector[neighbors.first] - var1_vector[neighbors.first - 1]);

            const auto eta = (val2 - var2_vector[neighbors.second - 1]) /
                             (var2_vector[neighbors.second] - var2_vector[neighbors.second - 1]);

            // use these coordinates for a bilinear interpolation
            const double interpolated_value =
                    (1 - xi) * (1 - eta) * eos_data[neighbors.first - 1][neighbors.second - 1][property_index] +
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

// a class to handle the specific needs of energy interpolation based on energy and entropy
// this class assumes that density is gridded and related to energy.  entropy is not gridded and therefore not directly related
class RestrictedBilinearInterpolation {
private:
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

    // basic function for calculating directional 1D distances
    static double calc_distance(double &given_value, double &sample_value) {
        double distance = given_value - sample_value;
        return distance;
    }

    // function that returns the bounds of the repeated densities that compose the nearest neighbors
    // specify "upper" or "lower" bounds for the "upper" or "lower" neighbors
    // these indices will be used to further restrict the var2/var3 vectors for nearest-neighbor searches
    static std::pair<unsigned int, unsigned int> get_neighbors_restricted_indices(
            const std::vector<double> &v,
            double value,
            unsigned int grid_length,
            const std::string &bound
    ) {
        unsigned int b1 = 0; // lower bound index
        unsigned int b2 = 0; // upper bound index
        int v_size = v.size();

        if (bound == "upper") {
            for (unsigned int a = 0; a < v.size(); a = a + 1) {
                double val = v[a];
                if (val >= value) {
                    b1 = a;
                    b2 = a + grid_length;
                    if (b2 > v.size() - 1) {
                        b2 = v.size() - 1;
                    }
                    break;
                }
            }
        } else {
            std::vector<double> reversed_vector(
                    v.size()); // declare a new vector which will be the reversed original vector given in the function
            std::reverse_copy(std::begin(v), std::end(v), std::begin(reversed_vector)); // reverse the vector

            for (unsigned int a = 0; a < reversed_vector.size(); a = a + 1) {
                double val = reversed_vector[a];
                if (val <= value) {
//                    b2 = v_size - (a - 1) - grid_length - 1;
//                    b1 = v_size - (a - 1) - (2 * grid_length);
                    b1 = (v_size - 1) - a - (grid_length - 1);
                    b2 = b1 + grid_length;
                    if (b1 < 0) {
                        b1 = 0;
                    }
                    if (b2 == 0) {
                        b2 = grid_length;
                    }
                    break;
                }
            }
        }

        if (b1 == 0 && b2 == 0) {
            b1 = 0;
            b2 = grid_length;
        }


        return std::make_pair(b1, b2);
    };

    // function that returns the 4 nearest var2 neighbors
    // these neighbors are within the upper and lower var1 bounds
    static std::pair<int, int> get_var2_neighbors(
            std::pair<unsigned int, unsigned int> &restricted_indices,
            const std::vector<double> &v,
            double &given_var2,
            unsigned int &grid_length
    ) {
        std::vector<double> var2_v_restricted = slice(v, restricted_indices.first, restricted_indices.second);

        bool initial_calc = true;
        double min_distance = 0;
        unsigned int min_distance_index = 0;

        for (unsigned int a = 0; a < var2_v_restricted.size(); a = a + 1) {
            double val = var2_v_restricted[a];
            double distance = calc_distance(given_var2, val);
            if (initial_calc) {
                min_distance = distance;
                min_distance_index = a;
                initial_calc = false;
            } else if (abs(distance) < abs(min_distance)) {
                min_distance = distance;
                min_distance_index = a;
            }
        }
        if (min_distance < 0) {
            if (min_distance_index <= 0) {
                return std::make_pair(min_distance_index, min_distance_index + 1);
            } else if (min_distance_index + 1 == grid_length) {
                return std::make_pair(min_distance_index - 2, min_distance_index - 1);
            } else {
                return std::make_pair(min_distance_index - 1, min_distance_index);
            }
        } else if (min_distance > 0) {
            if (min_distance_index + 1 == grid_length) {
                return std::make_pair(min_distance_index - 2, min_distance_index - 1);
            } else {
                return std::make_pair(min_distance_index, min_distance_index + 1);
            }
        } else {
            if (min_distance_index >= grid_length - 1) {
                return std::make_pair(min_distance_index - 2, min_distance_index - 1);
            } else if (min_distance_index < 0) {
                return std::make_pair(min_distance_index + 1, min_distance_index + 2);
            } else {
                return std::make_pair(min_distance_index, min_distance_index + 1);
            }
        }
    };

    static double interpolate_restricted(
            double &p1, // lower var1 neighbor value
            double &p2, // upper var1 neighbor value
            double &s11, // lower var1 neighbor lower var2 neighbor value
            double &s12, // lower var1 neighbor upper var2 neighbor value
            double &s21, // upper var1 neighbor lower var2 neighbor value
            double &s22, // upper var1 neighbor upper var2 neighbor value
            double &u11, // lower var1 neighbor lower var3 neighbor value
            double &u12, // lower var1 neighbor upper var3 neighbor value
            double &u21, // upper var1 neighbor lower var3 neighbor value
            double &u22, // upper var1 neighbor upper var3 neighbor value
            double &val1, // the var1 value being interpolated
            double &val2 // the var2 value being interpolated
    ) {

        // perform a series of 3 linear interpolations in order to arrive at the final interpolated var3 value
        double u1 = 0.0;
        double u2 = 0.0;
        double u = 0.0;
        if (s11 == s12) {
            u1 = u11;
        } else {
            u1 = linear_interpolate(s11, s12, val2, u11, u12);
        }
        if (s21 == s22) {
            u2 = u21;
        } else {
            u2 = linear_interpolate(s21, s22, val2, u21, u22);
        }
        if (p1 == p2) {
            u = u1;
        } else {
            u = linear_interpolate(p1, p2, val1, u1, u2);
        }

        return u;
    };

    // a class for performing simple linear interpolations
    static double linear_interpolate(
            double &x1,
            double &x2,
            double &x,
            double &q1,
            double &q2
    ) {
        double f = (((x2 - x) / (x2 - x1)) * q1) + (((x - x1) / (x2 - x1)) * q2);

        return f;
    };

public:

    // the public interpolation function for forming the var3 interpolation
    static double interpolate(
            double val1,
            double val2,
            const std::vector<double> &var1_vector,
            const std::vector<double> &var2_vector,
            const std::vector<double> &val3_vector,
            const unsigned int property_index,
            const std::vector<std::vector<std::array<double, 6>>> &eos_data,
            unsigned int grid_length
    ) {

        std::pair<unsigned int, unsigned int> neighbors_lower;
        std::pair<unsigned int, unsigned int> neighbors_upper;
        std::pair<unsigned int, unsigned int> restricted_index_pair_lower_bound = get_neighbors_restricted_indices(
                var1_vector, val1, grid_length, "lower");
        std::pair<unsigned int, unsigned int> restricted_index_pair_upper_bound = get_neighbors_restricted_indices(
                var1_vector, val1, grid_length, "upper");
        neighbors_lower = get_var2_neighbors(restricted_index_pair_lower_bound, var2_vector, val2,
                                             grid_length);  // lower var2 neighbors
        neighbors_upper = get_var2_neighbors(restricted_index_pair_upper_bound, var2_vector, val2,
                                             grid_length); // upper var2 neighbors


        double var1_vector_lower_neighbor = slice(
                var1_vector, restricted_index_pair_lower_bound.first, restricted_index_pair_lower_bound.second)[0];
        double var1_vector_upper_neighbor = slice(
                var1_vector, restricted_index_pair_upper_bound.first, restricted_index_pair_upper_bound.second)[0];
        std::vector<double> var2_vector_restricted_lower = slice(var2_vector,
                                                                 restricted_index_pair_lower_bound.first,
                                                                 restricted_index_pair_lower_bound.second);
        std::vector<double> var2_vector_restricted_upper = slice(var2_vector,
                                                                 restricted_index_pair_upper_bound.first,
                                                                 restricted_index_pair_upper_bound.second);
        std::vector<double> val3_vector_restricted_lower = slice(val3_vector,
                                                                 restricted_index_pair_lower_bound.first,
                                                                 restricted_index_pair_lower_bound.second);
        std::vector<double> val3_vector_restricted_upper = slice(val3_vector,
                                                                 restricted_index_pair_upper_bound.first,
                                                                 restricted_index_pair_upper_bound.second);

//        std::cout << "********************" << std::endl;
//        std::cout << val1 << std::endl;
//        std::cout << val2 << std::endl;
//        std::cout << var1_vector_lower_neighbor << std::endl;
//        std::cout << var1_vector_upper_neighbor << std::endl;
//        std::cout << var2_vector_restricted_lower[neighbors_lower.first] << std::endl;
//        std::cout << var2_vector_restricted_lower[neighbors_lower.second] << std::endl;
//        std::cout << var2_vector_restricted_lower[neighbors_upper.first] << std::endl;
//        std::cout << var2_vector_restricted_lower[neighbors_upper.second] << std::endl;
//        std::cout << val3_vector_restricted_lower[neighbors_lower.first] << std::endl;
//        std::cout << val3_vector_restricted_lower[neighbors_lower.second] << std::endl;
//        std::cout << val3_vector_restricted_lower[neighbors_upper.first] << std::endl;
//        std::cout << val3_vector_restricted_lower[neighbors_upper.second] << std::endl;
//        std::cout << val1 << std::endl;
//        std::cout << val2 << std::endl;
//        std::cout << "~~~~~~~~~~~~~~~~~~~~" << std::endl;

        const double interpolated_value = interpolate_restricted(
                var1_vector_lower_neighbor,
                var1_vector_upper_neighbor,
                var2_vector_restricted_lower[neighbors_lower.first],
                var2_vector_restricted_lower[neighbors_lower.second],
                var2_vector_restricted_upper[neighbors_upper.first],
                var2_vector_restricted_upper[neighbors_upper.second],
                val3_vector_restricted_lower[neighbors_lower.first],
                val3_vector_restricted_lower[neighbors_lower.second],
                val3_vector_restricted_upper[neighbors_upper.first],
                val3_vector_restricted_upper[neighbors_upper.second],
                val1,
                val2
        );

        return interpolated_value;
    };
};
