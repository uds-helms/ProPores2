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
#include "vector.h"
#include <ostream>

class VectorTest : public ::testing::Test {
public:
    VectorTest() :
            int_zero(),
            int_one(1, 1, 1),
            int_one_two_three(1, 2, 3),
            double_zero(),
            double_one(1.0, 1.0, 1.0),
            double_one_two_three(1.5, 2.5, 3.5) {}

    const Vec<int> int_zero;
    const Vec<int> int_one;
    const Vec<int> int_one_two_three;
    const Vec<double> double_zero;
    const Vec<double> double_one;
    const Vec<double> double_one_two_three;
};

TEST_F(VectorTest, empty_constructor) {
    EXPECT_EQ(0, int_zero.x);
    EXPECT_EQ(0, int_zero.y);
    EXPECT_EQ(0, int_zero.z);
    EXPECT_EQ(0.0, double_zero.x);
    EXPECT_EQ(0.0, double_zero.y);
    EXPECT_EQ(0.0, double_zero.z);
}

TEST_F(VectorTest, unit_constructor) {
    EXPECT_EQ(1, int_one.x);
    EXPECT_EQ(1, int_one.y);
    EXPECT_EQ(1, int_one.z);
    EXPECT_EQ(1.0, double_one.x);
    EXPECT_EQ(1.0, double_one.y);
    EXPECT_EQ(1.0, double_one.z);
}

TEST_F(VectorTest, one_two_three_constructor) {
    EXPECT_EQ(1, int_one_two_three.x);
    EXPECT_EQ(2, int_one_two_three.y);
    EXPECT_EQ(3, int_one_two_three.z);
    EXPECT_EQ(1.5, double_one_two_three.x);
    EXPECT_EQ(2.5, double_one_two_three.y);
    EXPECT_EQ(3.5, double_one_two_three.z);
}

TEST_F(VectorTest, copy_construtor) {
    Vec<int> vec_int = Vec<int>(int_one_two_three);
    EXPECT_EQ(Vec<int>(1, 2, 3), vec_int);
    Vec<double> vec_double = Vec<double>(double_one_two_three);
    EXPECT_EQ(Vec<double>(1.5, 2.5, 3.5), vec_double);
}

TEST_F(VectorTest, equality) {
    Vec<int> int_vec = Vec<int>(1, 2, 3);
    EXPECT_EQ(int_one_two_three, int_vec);
    EXPECT_NE(int_one, int_vec);
    EXPECT_NE(int_one, int_zero);
    Vec<double> double_vec = Vec<double>(1.5, 2.5, 3.5);
    EXPECT_EQ(double_one_two_three, double_vec);
    EXPECT_NE(double_one, double_vec);
    EXPECT_NE(double_one, double_zero);
}

TEST_F(VectorTest, plus_operator) {
    EXPECT_EQ(int_one, int_one + int_zero);
    EXPECT_EQ(int_one, int_zero + 1);
    EXPECT_EQ(Vec<int>(2, 3, 4), int_one + int_one_two_three);
    EXPECT_EQ(double_one, double_one + double_zero);
    EXPECT_EQ(double_one, double_zero + 1);
    EXPECT_EQ(Vec<double>(2.5, 3.5, 4.5), double_one + double_one_two_three);
    Vec<int> vec = Vec<int>(1, 1, 1);
    vec += int_one;
    EXPECT_EQ(Vec<int>(2, 2, 2), vec);
    EXPECT_EQ(Vec<int>(1, 1, 1), int_one);
    vec += 5;
    EXPECT_EQ(Vec<int>(7, 7, 7), vec);
}

TEST_F(VectorTest, minus_operator) {
    EXPECT_EQ(int_one, int_one - int_zero);
    EXPECT_EQ(int_zero, int_one - 1);
    EXPECT_EQ(Vec<int>(0, 1, 2), int_one_two_three - int_one);
    EXPECT_EQ(double_one, double_one - double_zero);
    EXPECT_EQ(double_zero, double_one - 1);
    EXPECT_EQ(Vec<double>(0.5, 1.5, 2.5), double_one_two_three - double_one);
    Vec<int> vec = Vec<int>(2, 2, 2);
    vec -= int_one;
    EXPECT_EQ(Vec<int>(1, 1, 1), vec);
    EXPECT_EQ(Vec<int>(1, 1, 1), int_one);
    vec -= 5;
    EXPECT_EQ(Vec<int>(-4, -4, -4), vec);
}

TEST_F(VectorTest, multiply_operator) {
    EXPECT_EQ(int_zero, int_one * 0);
    EXPECT_EQ(int_zero, int_zero * 1);
    EXPECT_EQ(Vec<int>(2, 2, 2), int_one * 2);
    EXPECT_EQ(Vec<int>(2, 4, 6), int_one_two_three * 2);
    EXPECT_EQ(Vec<int>(-2, -4, -6), int_one_two_three * (-2));
    EXPECT_EQ(double_zero, double_one * 0);
    EXPECT_EQ(double_zero, double_zero * 1);
    EXPECT_EQ(Vec<double>(2.0, 2.0, 2.0), double_one * 2);
    EXPECT_EQ(Vec<double>(3, 5, 7), double_one_two_three * 2);
    EXPECT_EQ(Vec<double>(-3, -5, -7), double_one_two_three * (-2));
    Vec<int> vec = Vec<int>(1, 1, 1);
    vec *= 5;
    EXPECT_EQ(Vec<int>(5, 5, 5), vec);
}

TEST_F(VectorTest, division_operator) {
    EXPECT_ANY_THROW(int_one / 0);
    EXPECT_EQ(int_zero, int_zero / 1);
    EXPECT_EQ(int_zero, int_one / 2);
    EXPECT_EQ(Vec<int>(0, 1, 1), int_one_two_three / 2);
    EXPECT_EQ(Vec<int>(0, -1, -1), int_one_two_three / (-2));
    EXPECT_ANY_THROW(double_one / 0);
    EXPECT_EQ(double_zero, double_zero / 1);
    EXPECT_EQ(Vec<double>(0.5, 0.5, 0.5), double_one / 2);
    EXPECT_EQ(Vec<double>(0.75, 1.25, 1.75), double_one_two_three / 2);
    EXPECT_EQ(Vec<double>(-0.75, -1.25, -1.75), double_one_two_three / (-2));
    Vec<int> vec = Vec<int>(1, 1, 1);
    vec /= 5;
    EXPECT_EQ(Vec<int>(0, 0, 0), vec);
    EXPECT_ANY_THROW(vec /= 0);
    Vec<double> vec_double = Vec<double>(1.0, 1.0, 1.0);
    vec_double /= 5.0;
    EXPECT_EQ(Vec<double>(1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0), vec_double);
    EXPECT_ANY_THROW(vec_double /= 0.0);
}

TEST_F(VectorTest, squared_length) {
    EXPECT_EQ(0.0, int_zero.length_squared());
    EXPECT_EQ(3.0, int_one.length_squared());
    EXPECT_EQ(14.0, int_one_two_three.length_squared());
    EXPECT_EQ(0.0, double_zero.length_squared());
    EXPECT_EQ(3.0, double_one.length_squared());
    EXPECT_EQ(22.0, Vec<double>(3.0, 2.0, 3.0).length_squared());

}

TEST_F(VectorTest, length) {
    EXPECT_EQ(0.0, int_zero.length());
    EXPECT_EQ(sqrt(3.0), int_one.length());
    EXPECT_EQ(sqrt(14.0), int_one_two_three.length());
    EXPECT_EQ(0.0, double_zero.length());
    EXPECT_EQ(sqrt(3.0), double_one.length());
    EXPECT_EQ(sqrt(22.0), Vec<double>(3.0, 2.0, 3.0).length());
}

TEST_F(VectorTest, unit) {
    // prevent no-discard compiler warning
    Vec<double> b;
    EXPECT_ANY_THROW(b = int_zero.unit());
    EXPECT_EQ(Vec<double>(1 / sqrt(3.0), 1 / sqrt(3.0), 1 / sqrt(3)), int_one.unit());
    EXPECT_EQ(Vec<double>(1 / sqrt(14.0), 2 / sqrt(14.0), 3 / sqrt(14.0)), int_one_two_three.unit());
    EXPECT_ANY_THROW(b = double_zero.unit());
    EXPECT_EQ(Vec<double>(1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)), double_one.unit());
    EXPECT_EQ(Vec<double>(1.5 / sqrt(20.75), 2.5 / sqrt(20.75), 3.5 / sqrt(20.75)), double_one_two_three.unit());
}

TEST_F(VectorTest, str) {
    std::ostringstream oss;
    oss << double_one_two_three;
    EXPECT_EQ("(1.5, 2.5, 3.5)", oss.str());
}

TEST_F(VectorTest, distance_squared) {
    EXPECT_EQ(14.0, distance_squared(int_zero, int_one_two_three));
    EXPECT_EQ(5.0, distance_squared(int_one, int_one_two_three));
    EXPECT_EQ(22.0, distance_squared(Vec<double>(3.0, 2.0, 3.0), double_zero));
}

TEST_F(VectorTest, distance) {
    EXPECT_EQ(sqrt(14.0), distance(int_zero, int_one_two_three));
    EXPECT_EQ(sqrt(5.0), distance(int_one, int_one_two_three));
    EXPECT_EQ(sqrt(22.0), distance(Vec<double>(3.0, 2.0, 3.0), double_zero));
}

TEST_F(VectorTest, dot_product) {
    EXPECT_EQ(0.0, dot_product(int_zero, int_one_two_three));
    EXPECT_EQ(6.0, dot_product(int_one, int_one_two_three));
    EXPECT_EQ(14.0, dot_product(int_one_two_three, int_one_two_three));
    EXPECT_EQ(7.5, dot_product(double_one, double_one_two_three));
    EXPECT_EQ(0.0, dot_product(double_zero, double_one_two_three));
}

TEST_F(VectorTest, cross_product) {
    EXPECT_EQ(Vec<int>(0, 0, 0), cross_product(int_zero, int_one_two_three));
    EXPECT_EQ(Vec<int>(1, -2, 1), cross_product(int_one, int_one_two_three));
    EXPECT_EQ(Vec<int>(0, 0, 0), cross_product(int_one_two_three, int_one_two_three));
    EXPECT_EQ(Vec<double>(1.0, -2.0, 1.0), cross_product(double_one, double_one_two_three));
    EXPECT_EQ(Vec<double>(0.0, 0.0, 0.0), cross_product(double_zero, double_one_two_three));
}

TEST_F(VectorTest, cosine_squared) {
    EXPECT_ANY_THROW(cosine_squared(int_zero, int_one_two_three));
    EXPECT_ANY_THROW(cosine_squared(int_one_two_three, int_zero));
    EXPECT_ANY_THROW(cosine_squared(double_zero, double_one_two_three));
    EXPECT_ANY_THROW(cosine_squared(double_one_two_three, double_zero));
    EXPECT_EQ(1.0, cosine_squared(int_one, int_one));
    EXPECT_EQ(1.0, cosine_squared(double_one, double_one));
    EXPECT_EQ(pow(6, 2) / (3 * 14), cosine_squared(int_one, int_one_two_three));
}

TEST_F(VectorTest, access) {
    EXPECT_ANY_THROW(int_one_two_three[3]);
    EXPECT_ANY_THROW(int_one_two_three[-1]);
    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(0, int_zero[i]);
        EXPECT_EQ(1, int_one[i]);
        EXPECT_EQ(i + 1, int_one_two_three[i]);
        EXPECT_EQ(0.0, double_zero[i]);
        EXPECT_EQ(1.0, double_one[i]);
        EXPECT_EQ(i + 1.5, double_one_two_three[i]);
    }
}

TEST(vector_tests, almost_perpendicular) {
    Vec<int> v1(Vec<int>(10, 0, 0) - Vec<int>(0, 0, 0));
    Vec<int> v2(Vec<int>(5, 5, 0) - Vec<int>(5, -5, 0));
    EXPECT_TRUE(almost_perpendicular(v1, v2));
    Vec<int> v3(Vec<int>(5, 5, 5) - Vec<int>(5, -5, 5));
    EXPECT_TRUE(almost_perpendicular(v1, v3));
    EXPECT_FALSE(almost_perpendicular(v2, v3));
}

TEST(vector_tests, closest) {
    Vec<double> v1(Vec<double>(10, 0, 0) - Vec<double>(0, 0, 0));
    Vec<double> v2(Vec<double>(5, 5, 0) - Vec<double>(5, -5, 0));
    Vec<double> v3(Vec<double>(5, 5, 5) - Vec<double>(5, -5, 5));
    auto res = closest(Vec<double>(10, 0, 0), Vec<double>(0, 0, 0), Vec<double>(5, 5, 0), Vec<double>(5, -5, 0));
    Vec<double> p1 = std::get<0>(res);
    Vec<double> p2 = std::get<1>(res);
    double d = std::get<2>(res);
    EXPECT_NEAR(0, d, 0.000001);
    for (const Vec<double> &v: {p1, p2}) {
        EXPECT_NEAR(5, v.x, 0.000001);
        EXPECT_NEAR(0, v.y, 0.000001);
        EXPECT_NEAR(0, v.z, 0.000001);
    }
    res = closest(Vec<double>(5, 5, 5), Vec<double>(5, -5, 5), Vec<double>(10, 0, 0), Vec<double>(0, 0, 0));
    p1 = std::get<0>(res);
    p2 = std::get<1>(res);
    d = std::get<2>(res);
    EXPECT_NEAR(5, d, 0.000001);
    EXPECT_NEAR(5, p1.x, 0.000001);
    EXPECT_NEAR(0, p1.y, 0.000001);
    EXPECT_NEAR(5, p1.z, 0.000001);
    EXPECT_NEAR(5, p2.x, 0.000001);
    EXPECT_NEAR(0, p2.y, 0.000001);
    EXPECT_NEAR(0, p2.z, 0.000001);
}

TEST(vector_tests, vec_in_interval) {
    EXPECT_TRUE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(0, 0, 0), Vec<int>(0, 0, 0)));
    EXPECT_FALSE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(1, 0, 0), Vec<int>(0, 0, 0)));
    EXPECT_FALSE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(0, 1, 0), Vec<int>(0, 0, 0)));
    EXPECT_FALSE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(0, 0, 1), Vec<int>(0, 0, 0)));

    EXPECT_TRUE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(0, 0, 0), Vec<int>(1, 1, 1)));
    EXPECT_TRUE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(1, 0, 0), Vec<int>(1, 1, 1)));
    EXPECT_TRUE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(0, 1, 0), Vec<int>(1, 1, 1)));
    EXPECT_TRUE(vec_in_interval(Vec<int>(0, 0, 0), Vec<int>(0, 0, 1), Vec<int>(1, 1, 1)));
}
