/*
 * Copyright (C) 2020 Markus Hollander (markus.hollander@bioinformatik.uni-saarland.de). All rights reserved.
 *
 * This file is part of PROPORES 2.0. PROPORES was originally developed by Po-Hsien Lee (2011) in Perl under the terms
 * of the GNU General Public License version 3 or higher. The approach was published as:
 * Lee, PH, Helms, V (2012). Identifying continuous pores in protein structures with PROPORES by computational
 * repositioning of gating residues. Proteins, 80, 2:421-32. https://www.ncbi.nlm.nih.gov/pubmed/22095919
 *
 * PROPORES 2.0 is a C++ implementation of the approach outlined in Lee (2012) by Markus Hollander and Moomal Aziz.
 * It is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the license, or (at your option) any later version.
 *
 * PROPORES 2.0 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org./licenses/>.
*/

#ifndef PROPORES_VECTOR_H
#define PROPORES_VECTOR_H

#include <cmath>
#include <tuple>
#include <ostream>
#include <exception>
#include <algorithm>

template<typename T>
struct Vec {
    T x;
    T y;
    T z;

    Vec() : x(), y(), z() {};

    Vec(T const x_in, T const y_in, T const z_in) : x(x_in), y(y_in), z(z_in) {};

    Vec(const Vec<T> &vec) : x(vec.x), y(vec.y), z(vec.z) {};

    T operator[](int const i) const {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                throw std::runtime_error("Vec: index must be 0 <= i <= 2.");
        }
    }

    // vector addition
    Vec<T> operator+(const Vec<T> &vec) const { return Vec<T>(x + vec.x, y + vec.y, z + vec.z); }

    // scalar addition
    Vec<T> operator+(T const s) const { return Vec<T>(x + s, y + s, z + s); }

    // vector addition assignment
    Vec<T> operator+=(const Vec<T> &vec) {
        x += vec.x;
        y += vec.y;
        z += vec.z;
        return *this;
    }

    // scalar addition assignment
    Vec<T> operator+=(const T s) {
        x += s;
        y += s;
        z += s;
        return *this;
    }

    // vector subtraction
    Vec<T> operator-(const Vec<T> &vec) const { return Vec<T>(x - vec.x, y - vec.y, z - vec.z); }

    // scalar subtraction
    Vec<T> operator-(T const s) const { return Vec<T>(x - s, y - s, z - s); }

    // vector subtraction assignment
    Vec<T> operator-=(const Vec<T> &vec) {
        x -= vec.x;
        y -= vec.y;
        z -= vec.z;
        return *this;
    }

    // scalar subtraction assignment
    Vec<T> operator-=(const T s) {
        x -= s;
        y -= s;
        z -= s;
        return *this;
    }

    // scalar multiplication
    Vec<T> operator*(T const s) const { return Vec<T>(x * s, y * s, z * s); }

    // scalar multiplication assignment
    Vec<T> operator*=(const T s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    // scalar division
    Vec<T> operator/(T const s) const {
        if (s == 0) throw std::invalid_argument("Vector division by zero-scalar.");
        return Vec<T>(x / s, y / s, z / s);
    }

    // scalar division assignment
    Vec<T> operator/=(const T s) {
        if (s == 0) throw std::invalid_argument("Vector division by zero-scalar.");
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    // equality
    bool operator==(const Vec<T> &vec) const { return x == vec.x && y == vec.y && z == vec.z; }

    bool operator!=(const Vec<T> &vec) const { return x != vec.x || y != vec.y || z != vec.z; }

    // length
    [[nodiscard]] double length() const { return sqrt(length_squared()); }

    [[nodiscard]] T length_squared() const { return T(pow(x, 2) + pow(y, 2) + pow(z, 2)); }

    // transforms the vector into a directional vector of length 1
    [[nodiscard]] Vec<double> unit() const {
        double norm = length();
        if (norm == 0.0) throw std::runtime_error("Null-vector cannot be converted to unit-vector.");
        return Vec<double>(x / norm, y / norm, z / norm);
    }
};

template<typename T>
std::ostream &operator<<(std::ostream &os, const Vec<T> &vec) {
    return os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";  // whatever needed to print vec to os
}

template<typename A, typename B>
double distance_squared(const Vec<A> &a, const Vec<B> &b) {
    return pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
}

template<typename A, typename B>
double distance(const Vec<A> &a, const Vec<B> &b) { return sqrt(distance_squared(a, b)); }

template<typename A, typename B>
double dot_product(const Vec<A> &a, const Vec<B> &b) { return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }

template<typename A, typename B>
double cosine_squared(const Vec<A> &a, const Vec<B> &b) {
    // cosine similarity: (a * b) / (a.length() * b.length())
    if (a.length_squared() == 0.0) throw std::runtime_error("Vector A has length 0.");
    if (b.length_squared() == 0.0) throw std::runtime_error("Vector B has length 0.");
    return pow(dot_product(a, b), 2) / (a.length_squared() * b.length_squared());
}

// cos(85)**2 = cos(95)**2 = 0.007596124
// therefore cos(x) <= 0.007596124 means it's almost perpendicular
template<typename A, typename B>
bool almost_perpendicular(const Vec<A> &a, const Vec<B> &b) { return cosine_squared(a, b) <= 0.007596124; }

template<typename T>
// compute the shortest distance between two 3D lines (point_1 -> point_2 and point_3 -> point_4)
// by identifying point_a on the first line and point_b on the second line with minimal distance
// from: http://paulbourke.net/geometry/pointlineplane/lineline.c
std::tuple<Vec<double>, Vec<double>, double> closest(const Vec<T> &point_1, const Vec<T> &point_2,
                                                     const Vec<T> &point_3, const Vec<T> &point_4) {
    Vec<T> vec_13 = point_1 - point_3;
    Vec<T> vec_43 = point_4 - point_3;
    Vec<T> vec_21 = point_2 - point_1;
    // compute the dot products
    double dot_13_43 = dot_product(vec_13, vec_43);
    double dot_43_21 = dot_product(vec_43, vec_21);
    double dot_13_21 = dot_product(vec_13, vec_21);
    double dot_43 = dot_product(vec_43, vec_43);
    double dot_21 = dot_product(vec_21, vec_21);

    double mu_a = (dot_13_43 * dot_43_21 - dot_13_21 * dot_43) / (dot_21 * dot_43 - dot_43_21 * dot_43_21);
    double mu_b = (dot_13_43 + dot_43_21 * mu_a) / dot_43;

    Vec<double> point_a = point_1 + vec_21 * mu_a;
    Vec<double> point_b = point_3 + vec_43 * mu_b;
    double dist = distance(point_a, point_b);

    return std::make_tuple(point_a, point_b, dist);
}

template<typename T>
bool vec_in_interval(const Vec<T> &a, const Vec<T> &vec, const Vec<T> &b) {
    if (vec.x < std::min(a.x, b.x) || vec.x > std::max(a.x, b.x)) return false;
    if (vec.y < std::min(a.y, b.y) || vec.y > std::max(a.y, b.y)) return false;
    return vec.z >= std::min(a.z, b.z) && vec.z <= std::max(a.z, b.z);
}

template<typename T>
Vec<T> cross_product(const Vec<T> &a, const Vec<T> &b) {
    return Vec<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

#endif //PROPORES_VECTOR_H
