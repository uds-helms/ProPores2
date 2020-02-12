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

#include "gtest/gtest.h"
#include "angles.h"
#include <iostream>

TEST(angles, dihedral) {
    Vec<double> a(-10, 0, 0);
    Vec<double> b(0, 0, 0);
    Vec<double> c(0, 0, -10);
    Vec<double> d(0, 10, -10);

    double res_1 = dihedral_angle(a, b, c, d);
    double res_2 = dihedral_angle_2(a, b, c, d);
    EXPECT_NEAR(res_1, res_2, 0.0001);
    EXPECT_NEAR(radian_to_degree(res_1), 90, 0.0001);

    d = Vec<double>(0, -10, -10);
    res_1 = dihedral_angle(a, b, c, d);
    res_2 = dihedral_angle_2(a, b, c, d);
    EXPECT_NEAR(res_1, res_2, 0.0001);
    EXPECT_NEAR(radian_to_degree(res_1), 270, 0.0001);
}

TEST(angles, quaternion_rotation_matrix) {
    double rotation_angle = 30;
    Vec<double> rotation_axis(1, 1, 1);
    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, rotation_axis);
    EXPECT_NEAR(0.9107, row_1.x, 0.001);
    EXPECT_NEAR(0.9107, row_2.y, 0.001);
    EXPECT_NEAR(0.9107, row_3.z, 0.001);
    EXPECT_NEAR(0.3333, row_1.z, 0.001);
    EXPECT_NEAR(0.3333, row_2.x, 0.001);
    EXPECT_NEAR(0.3333, row_3.y, 0.001);
    EXPECT_NEAR(-0.244, row_1.y, 0.001);
    EXPECT_NEAR(-0.244, row_2.z, 0.001);
    EXPECT_NEAR(-0.244, row_3.x, 0.001);

    Vec<double> v(2, 1, 3);
    Vec<double> v_rotated(dot_product(row_1, v), dot_product(row_2, v), dot_product(row_3, v));
    EXPECT_NEAR(2.5773, v_rotated.x, 0.001);
    EXPECT_NEAR(0.8453, v_rotated.y, 0.001);
    EXPECT_NEAR(2.5774, v_rotated.z, 0.001);
}

TEST(angles, quaternion_rotation_matrix_2) {
    double rotation_angle = 30;
    Vec<double> rotation_axis(1, 1, 1);
    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, rotation_axis);
    EXPECT_NEAR(0.9107, row_1.x, 0.001);
    EXPECT_NEAR(0.9107, row_2.y, 0.001);
    EXPECT_NEAR(0.9107, row_3.z, 0.001);
    EXPECT_NEAR(0.3333, row_1.z, 0.001);
    EXPECT_NEAR(0.3333, row_2.x, 0.001);
    EXPECT_NEAR(0.3333, row_3.y, 0.001);
    EXPECT_NEAR(-0.244, row_1.y, 0.001);
    EXPECT_NEAR(-0.244, row_2.z, 0.001);
    EXPECT_NEAR(-0.244, row_3.x, 0.001);

    Vec<double> v(2, 1, 3);
    Vec<double> v_rotated(dot_product(row_1, v), dot_product(row_2, v), dot_product(row_3, v));
    EXPECT_NEAR(2.5773, v_rotated.x, 0.001);
    EXPECT_NEAR(0.8453, v_rotated.y, 0.001);
    EXPECT_NEAR(2.5774, v_rotated.z, 0.001);
}

TEST(angles, quaternion_rotation_matrix_negative) {
    double rotation_angle = -30;
    Vec<double> rotation_axis(1, 1, 1);
    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, rotation_axis);
    EXPECT_NEAR(0.9107, row_1.x, 0.001);
    EXPECT_NEAR(0.9107, row_2.y, 0.001);
    EXPECT_NEAR(0.9107, row_3.z, 0.001);
    EXPECT_NEAR(-0.244, row_1.z, 0.001);
    EXPECT_NEAR(-0.244, row_2.x, 0.001);
    EXPECT_NEAR(-0.244, row_3.y, 0.001);
    EXPECT_NEAR(0.3333, row_1.y, 0.001);
    EXPECT_NEAR(0.3333, row_2.z, 0.001);
    EXPECT_NEAR(0.3333, row_3.x, 0.001);

    Vec<double> v(2, 1, 3);
    Vec<double> v_rotated(dot_product(row_1, v), dot_product(row_2, v), dot_product(row_3, v));
    EXPECT_NEAR(1.4227, v_rotated.x, 0.001);
    EXPECT_NEAR(1.4226, v_rotated.y, 0.001);
    EXPECT_NEAR(3.1547, v_rotated.z, 0.001);
}

TEST(angles, quaternion_rotation_matrix_negative_2) {
    double rotation_angle = -390;
    Vec<double> rotation_axis(1, 1, 1);
    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, rotation_axis);
    EXPECT_NEAR(0.9107, row_1.x, 0.001);
    EXPECT_NEAR(0.9107, row_2.y, 0.001);
    EXPECT_NEAR(0.9107, row_3.z, 0.001);
    EXPECT_NEAR(-0.244, row_1.z, 0.001);
    EXPECT_NEAR(-0.244, row_2.x, 0.001);
    EXPECT_NEAR(-0.244, row_3.y, 0.001);
    EXPECT_NEAR(0.3333, row_1.y, 0.001);
    EXPECT_NEAR(0.3333, row_2.z, 0.001);
    EXPECT_NEAR(0.3333, row_3.x, 0.001);

    Vec<double> v(2, 1, 3);
    Vec<double> v_rotated(dot_product(row_1, v), dot_product(row_2, v), dot_product(row_3, v));
    EXPECT_NEAR(1.4227, v_rotated.x, 0.001);
    EXPECT_NEAR(1.4226, v_rotated.y, 0.001);
    EXPECT_NEAR(3.1547, v_rotated.z, 0.001);
}

TEST(angles, integration) {
    Vec<double> a(-10, 0, 0);
    Vec<double> b(0, 0, 0);
    Vec<double> c(0, 0, -10);
    Vec<double> d(0, 10, -10);

    double res_1 = dihedral_angle(a, b, c, d);
    EXPECT_NEAR(radian_to_degree(res_1), 90, 0.0001);

    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(90, c - b);
    Vec<double> rotated_d(dot_product(row_1, d), dot_product(row_2, d), dot_product(row_3, d));
    EXPECT_NEAR(10, rotated_d.x, 0.001);
    EXPECT_NEAR(0, rotated_d.y, 0.001);
    EXPECT_NEAR(-10, rotated_d.z, 0.001);
    EXPECT_NEAR(180, abs(radian_to_degree(dihedral_angle(a, b, c, rotated_d))), 0.001);
}

TEST(angles, integration_2) {
    Vec<double> a(-10, 1, 1);
    Vec<double> b(0, 1, 1);
    Vec<double> c(0, 1, -10);
    Vec<double> d(0, 10, -10);
    double res_1 = dihedral_angle(a, b, c, d);
    EXPECT_NEAR(radian_to_degree(res_1), 90, 0.0001);
    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(90, c - b);
    Vec<double> tmp = d - c;
    Vec<double> d_rot = Vec<double>(dot_product(row_1, tmp), dot_product(row_2, tmp), dot_product(row_3, tmp)) + c;
    EXPECT_NEAR(9, d_rot.x, 0.001);
    EXPECT_NEAR(1, d_rot.y, 0.001);
    EXPECT_NEAR(-10, d_rot.z, 0.001);
    EXPECT_NEAR(180, abs(radian_to_degree(dihedral_angle(a, b, c, d_rot))), 0.001);
}

TEST(angles, integration_3) {
    Vec<double> a(0.472, 0.98, 16.392);
    Vec<double> b(1.206, 1.955, 15.591);
    Vec<double> c(2.676, 1.549, 15.454);
    Vec<double> d(2.909, 0.318, 14.598);
    double rotation_angle = -69.2557;

    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, c - b);
    Vec<double> tmp = d - c;
    Vec<double> d_rot = Vec<double>(dot_product(row_1, tmp), dot_product(row_2, tmp), dot_product(row_3, tmp)) + c;
    double angle = radian_to_degree(dihedral_angle(a, b, c, d_rot));
    if (angle > 300) {
        EXPECT_NEAR(360, angle, 0.001);
    } else {
        EXPECT_NEAR(0, angle, 0.001);
    }

}

TEST(angles, integration_4) {
    /*
     * MET: N=(0.472, 0.98, 16.392) CA=(1.206, 1.955, 15.591) CB=(2.676, 1.549, 15.454) CG=(2.909, 0.318, 14.598) SD=(4.644, 0.
        144, 14.121) CE=(4.7, 1.204, 12.68)
        0: old angle=69.2557 | rotation=200.744 | desired=270 | new angle=270
        1: old angle=165.051 | rotation=-45.0505 | desired=120 | new angle=219.723
        2: old angle=275.725 | rotation=-245.725 | desired=30 | new angle=296.418
     */
    Vec<double> a(0.472, 0.98, 16.392);
    Vec<double> b(1.206, 1.955, 15.591);
    Vec<double> c(2.676, 1.549, 15.454);
    Vec<double> d(2.909, 0.318, 14.598);
    Vec<double> e(4.644, 0.144, 14.121);
    Vec<double> f(4.7, 1.204, 12.68);
    std::vector<Vec<double>> vecs = {a, b, c, d, e, f};

    double angle_1 = radian_to_degree(dihedral_angle(a, b, c, d));
    EXPECT_NEAR(69.2557, angle_1, 0.001);
    double angle_2 = radian_to_degree(dihedral_angle(b, c, d, e));
    EXPECT_NEAR(165.051, angle_2, 0.001);
    double angle_3 = radian_to_degree(dihedral_angle(c, d, e, f));
    EXPECT_NEAR(275.725, angle_3, 0.001);

    double rotation_angle = 200.744;
    auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, c - b);
    Vec<double> tmp = d - c;
    d = Vec<double>(dot_product(row_1, tmp), dot_product(row_2, tmp), dot_product(row_3, tmp)) + c;
    tmp = e - c;
    e = Vec<double>(dot_product(row_1, tmp), dot_product(row_2, tmp), dot_product(row_3, tmp)) + c;
    tmp = f - c;
    f = Vec<double>(dot_product(row_1, tmp), dot_product(row_2, tmp), dot_product(row_3, tmp)) + c;
    EXPECT_NEAR(270, radian_to_degree(dihedral_angle(a, b, c, d)), 0.001);
    EXPECT_NEAR(angle_2, radian_to_degree(dihedral_angle(b, c, d, e)), 0.001);
    EXPECT_NEAR(angle_3, radian_to_degree(dihedral_angle(c, d, e, f)), 0.001);

    rotation_angle = -45.0505;
    auto[row_12, row_22, row_32] = quaternion_rotation_matrix(rotation_angle, d - c);
    tmp = e - d;
    e = Vec<double>(dot_product(row_12, tmp), dot_product(row_22, tmp), dot_product(row_32, tmp)) + d;
    tmp = f - d;
    f = Vec<double>(dot_product(row_12, tmp), dot_product(row_22, tmp), dot_product(row_32, tmp)) + d;
    EXPECT_NEAR(270, radian_to_degree(dihedral_angle(a, b, c, d)), 0.001);
    EXPECT_NEAR(120, radian_to_degree(dihedral_angle(b, c, d, e)), 0.001);
    EXPECT_NEAR(angle_3, radian_to_degree(dihedral_angle(c, d, e, f)), 0.001);

    rotation_angle = -245.725;
    auto[row_13, row_23, row_33] = quaternion_rotation_matrix(rotation_angle, e - d);
    tmp = f - e;
    f = Vec<double>(dot_product(row_13, tmp), dot_product(row_23, tmp), dot_product(row_33, tmp)) + e;
    EXPECT_NEAR(270, radian_to_degree(dihedral_angle(a, b, c, d)), 0.001);
    EXPECT_NEAR(120, radian_to_degree(dihedral_angle(b, c, d, e)), 0.001);
    EXPECT_NEAR(30, radian_to_degree(dihedral_angle(c, d, e, f)), 0.001);
}