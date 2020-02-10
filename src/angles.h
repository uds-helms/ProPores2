#ifndef PROPORES_ANGLES_H
#define PROPORES_ANGLES_H

#include <iostream>
#include "vector.h"

const double PI = 3.14159265359;

// convert an angle from degree to radian
double degree_to_radian(double angle);

// convert an angle from radian to degree, convert negative to positive degrees
double radian_to_degree(double angle);

// dihedral angle between the plane of points a and b, and the plane of points c and d
// both planes are connected by the vector between points b and c
// sightly faster due to fewer cross-products
template<typename T>
double dihedral_angle(const Vec<T> &a, const Vec<T> &b, const Vec<T> &c, const Vec<T> &d) {
    // vectors from b to a, b to c, and c to d
    Vec<T> ba = (b - a) * -1;
    Vec<T> bc_unit = (c - b).unit();
    Vec<T> cd = d - c;
    // ba and cd projected onto the plane perpendicular to bc_unit, respectively
    Vec<T> v = ba - bc_unit * dot_product(ba, bc_unit);
    Vec<T> w = cd - bc_unit * dot_product(cd, bc_unit);
    // the dihedral angle is the angle between v and w in the plane, which projected into Euclidean space is the
    // angle between the point (x, y) and the positive x-axis
    double x = dot_product(v, w);
    double y = dot_product(cross_product(bc_unit, v), w);
    // the resulting angle is in radian: -pi < angle <= pi
    return atan2(y, x);
}

// dihedral angle between the plane of points a and b, and the plane of points c and d
// both planes are connected by the vector between points b and c
// slightly slower due to more cross-products
template<typename T>
double dihedral_angle_2(const Vec<T> &a, const Vec<T> &b, const Vec<T> &c, const Vec<T> &d) {
    // vectors from b to a, b to c, and c to d
    Vec<T> ba = (b - a) * -1;
    Vec<T> bc = c - b;
    Vec<T> cd = d - c;
    // vectors normal to the first and second plane, as well as normal to the plane formed by the two normal vectors
    Vec<T> normal_1 = cross_product(ba, bc);
    Vec<T> normal_2 = cross_product(cd, bc);
    Vec<T> normal_3 = cross_product(normal_1, normal_2);
    // the dihedral angle is the angle between the normal vectors, which projected into Euclidean space is the
    // angle between the point (x, y) and the positive x-axis
    double x = dot_product(normal_1, normal_2);
    double y = dot_product(normal_3, bc) / bc.length();
    // the resulting angle is in radian: -pi < angle <= pi
    return atan2(y, x);
}

template<typename T>
std::tuple<Vec<double>, Vec<double>, Vec<double>> quaternion_rotation_matrix(const double angle_degree,
                                                                             const Vec<T> &rotation_axis) {
    // convert the rotation angle to radian and divide it by 2 since the following operations work with theta / 2
    double theta = degree_to_radian(angle_degree) / 2;
    // scalar part of the quaternion that describes the rotation
    double q_r = cos(theta);
    // vector part of the quaternion that describes the rotation axis
    Vec<double> q = rotation_axis.unit() * sin(theta);
    // normalisation factor
    double s = 2 / (pow(q_r, 2) + pow(q.x, 2) + pow(q.y, 2) + pow(q.z, 2));
    // rotation matrix cell components
    double xr = q.x * q_r;
    double yr = q.y * q_r;
    double zr = q.z * q_r;
    double xx = q.x * q.x;
    double yy = q.y * q.y;
    double zz = q.z * q.z;
    double xy = q.x * q.y;
    double xz = q.x * q.z;
    double yz = q.y * q.z;
    // construct the rotation matrix from the components
    Vec<double> row_1(1 - s * (yy + zz), s * (xy - zr), s * (xz + yr));
    Vec<double> row_2(s * (xy + zr), 1 - s * (xx + zz), s * (yz - xr));
    Vec<double> row_3(s * (xz - yr), s * (yz + xr), 1 - s * (xx + yy));
    // return the rows
    return std::make_tuple(row_1, row_2, row_3);
}

#endif //PROPORES_ANGLES_H
