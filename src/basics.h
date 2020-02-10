#ifndef PROPORES_BASICS_H
#define PROPORES_BASICS_H

#include <map>
#include <string>
#include <chrono>
#include <memory>
#include <vector>
#include <fstream>
#include <sstream>
#include <exception>
#include <algorithm>
#include <unordered_set>
#include "vector.h"

// construct an indent (whitespace) for the given indent level
std::string indent(size_t level);

// print a message
void print(const std::string &msg);

// print an indent for the specified level followed by the message
void print(size_t level, const std::string &msg);

// print an indent of the specified level followed by the message and the passed time in seconds
void print(size_t level, const std::chrono::time_point<std::chrono::high_resolution_clock> &start_time,
           const std::string &msg);

// remove all spaces from a given string
std::string remove_spaces(std::string str);

// strip whitespace from the right side of a string
std::string r_strip(std::string str);

// extract the substring up to the given length from the input string
std::string substr(const std::string &input_str, size_t start, size_t length);

// extract the substring up to the given length from the input string, and remove asterisks
std::string pdb_substr(const std::string &input_str, size_t start, size_t length);

double pdb_substr(const std::string &input_str, size_t start, size_t length, double default_value);

// split the string at whitespace
std::vector<std::string> split(const std::string &str);

// split the string at the given separator
std::vector<std::string> split(const std::string &str, char sep);

// split the words of the input string into rows such that each row is at most <width> characters long
std::vector<std::string> wrap(const std::string &str, size_t width);

// search box relative to a grid box of interest
struct RelativeSearchBox {
    // relative coordinates on the grid => can be negative
    Vec<int> grid;
    // squared distance from the search centre
    int length_sq;

    RelativeSearchBox(const int x, const int y, const int z) :
            grid(x, y, z),
            length_sq(grid.length_squared()) {}
};

// cube of relative search boxes that can be used to search around a grid box
struct RelativeSearchSpaces {
    std::map<size_t, std::vector<RelativeSearchBox>> space_map;

    // radius = number of relative grid boxes on each side
    // generate the search space if it does not exist yet, otherwise just return it
    std::vector<RelativeSearchBox> &operator[](const int radius) {
        // relative search space was already constructed
        if (space_map.count(radius)) return space_map.at(radius);
        // initialise new relative search space
        std::vector<RelativeSearchBox> search_space;
        // construct a relative cube with the given radius
        for (int x = -radius; x <= radius; x++) {
            for (int y = -radius; y <= radius; y++) {
                for (int z = -radius; z <= radius; z++) {
                    search_space.emplace_back(RelativeSearchBox(x, y, z));
                }
            }
        }
        // store the new relative search boxes and then return them
        space_map.insert(std::make_pair(radius, search_space));
        return space_map.at(radius);
    }
};

#endif //PROPORES_BASICS_H
