#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include "basics.h"

// construct an indent (whitespace) for the given indent level
std::string indent(const size_t level) { return std::string(level * 2, ' '); }

// print a message
void print(const std::string &msg) { std::cout << msg << std::endl; }

// print an indent for the specified level followed by the message
void print(size_t level, const std::string &msg) { std::cout << indent(level) << msg << std::endl; }

// print an indent of the specified level followed by the message and the passed time in seconds
void print(const size_t level, const std::chrono::time_point<std::chrono::high_resolution_clock> &start_time,
           const std::string &msg) {
    auto current_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
    std::cout << indent(level) << msg << " " << std::setprecision(2) << std::fixed << duration << "s" << std::endl;
}

// remove all spaces from a given string
std::string remove_spaces(std::string str) {
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    return str;
}

// strip whitespace from the right side of a string
std::string r_strip(std::string str) {
    str.erase(std::find_if(str.rbegin(), str.rend(), [](char c) { return !std::isspace(c); }).base(), str.end());
    return str;
}

// extract the substring up to the given length from the input string
std::string substr(const std::string &str, const size_t start, const size_t length) {
    // empty substring
    if (length < 1) return "";
    // the start of the substring is outside the string
    if (str.size() <= start) return "";
    // get the substring of the given length, or the rest of the string
    return str.substr(start, std::min(length, str.size() - start));
}

// extract the substring up to the given length from the input string, and remove asterisks
std::string pdb_substr(const std::string &input_str, const size_t start, const size_t length) {
    // get the substring and replace asterisks with spaces
    std::string str = substr(input_str, start, length);
    std::replace(str.begin(), str.end(), '*', ' ');
    // remove spaces before returning the modified string
    return remove_spaces(str);
}

double pdb_substr(const std::string &input_str, size_t start, size_t length, double default_value) {
    std::string str = pdb_substr(input_str, start, length);
    if (str.empty()) return default_value;
    return stod(str);
}

// split the string at whitespace
std::vector<std::string> split(const std::string &str) {
    std::vector<std::string> items;
    std::stringstream ss(str);
    std::string item;
    while (ss >> item) { items.push_back(item); }
    return items;
}

// split a string at the given separator
std::vector<std::string> split(const std::string &str, const char sep) {
    std::vector<std::string> items;
    std::stringstream ss(str);
    std::string item;
    while (getline(ss, item, sep)) { items.push_back(item); }
    return items;
}

// split the words of the input string into rows such that each row is at most <width> characters long
std::vector<std::string> wrap(const std::string &str, const size_t width) {
    std::vector<std::string> words = split(str);
    std::vector<std::string> rows;
    // construct the rows
    for (const std::string &word: words) {
        // create a new row for the first word and also for words that would push the current row past the width limit
        if (rows.empty() || rows[rows.size() - 1].length() + word.length() >= width) {
            rows.push_back(word);
            // otherwise just append the word to the current row
        } else {
            rows[rows.size() - 1].append(" " + word);
        }
    }
    return rows;
}
