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

#include <memory>
#include "grid.h"
#include "enums.h"
#include "reader.h"
#include "gtest/gtest.h"

class ProteinGridTest : public ::testing::Test {
public:
    ProteinGridTest() :
            small("test_files/simple_coord.pdb", 1.0, 1.0, false, false),
            small_05("test_files/simple_coord.pdb", 0.5, 1.0, false, false),
            kd("test_files/kd.pdb", 1.0, 1.0, false, false),
            g_per("test_files/perpendicular.pdb", 1.0, 1.0, false, false) {
        atom_1 = std::make_shared<Atom>(Atom{0, "ATOM      3  C   SER A  17    0.123456 1.123452.123456"});
        atom_2 = std::make_shared<Atom>(Atom{0, "ATOM      3  C   SER A  17    0.123456 1.123452.123456"});
        atom_3 = std::make_shared<Atom>(Atom{1, "ATOM      3  C   SER A  17    0.123456 1.123452.123456"});
    }

    ProteinGrid small;
    ProteinGrid small_05;
    ProteinGrid kd;
    ProteinGrid g_per;
    std::shared_ptr<Atom> atom_1;
    std::shared_ptr<Atom> atom_2;
    std::shared_ptr<Atom> atom_3;
};

TEST_F(ProteinGridTest, valid) {
    EXPECT_TRUE(small.valid(small.min));
    EXPECT_FALSE(small.valid(small.max));
    EXPECT_FALSE(small.valid(Vec<int>(6, 8, 6)));
    EXPECT_FALSE(small.valid(Vec<int>(7, 7, 6)));
    EXPECT_FALSE(small.valid(Vec<int>(6, 7, 7)));
    EXPECT_TRUE(small.valid(Vec<int>(-5, 4, 0)));
}

TEST_F(ProteinGridTest, box) {
    EXPECT_EQ(Vec<int>(3, 4, 5), small.box(Vec<double>(3.1, 4.1, 5.1)));
    EXPECT_EQ(Vec<int>(6, 8, 10), small_05.box(Vec<double>(3.1, 4.1, 5.1)));
}

TEST_F(ProteinGridTest, centre) {
    Vec<double> actual = small.centre(1, 1, 1);
    Vec<double> expected = Vec<double>(1.5, 1.5, 1.5);
    for (size_t i = 0; i < 3; i++) {
        EXPECT_NEAR(expected[i], actual[i], 0.0001);
    }

    actual = small.centre(0, 5, -1);
    expected = Vec<double>(0.5, 5.5, -0.5);
    for (size_t i = 0; i < 3; i++) {
        EXPECT_NEAR(expected[i], actual[i], 0.0001);
    }

    actual = small_05.centre(0, 5, -1);
    expected = Vec<double>(0.25, 2.75, -0.25);
    for (size_t i = 0; i < 3; i++) {
        EXPECT_NEAR(expected[i], actual[i], 0.0001);
    }
}

TEST_F(ProteinGridTest, index_and_at) {
    EXPECT_EQ(0, small.index(small.min));
    EXPECT_EQ(2562, small.index(small.max));

    EXPECT_EQ(0, small_05.index(small_05.min));
    EXPECT_EQ(18981, small_05.index(small_05.max));

    small.states = std::vector<BoxState>(27);
    small.height = 3;
    small.depth = 3;
    small.width = 3;
    small.boxes = 27;
    small.min = Vec<int>(-1, -1, -1);
    small.max = Vec<int>(2, 2, 2);

    EXPECT_EQ(0, small.index(-1, -1, -1));
    EXPECT_EQ(52, small.index(3, 3, 3));

    EXPECT_EQ(5, small.index(1, 0, -1));
    EXPECT_EQ(22, small.index(0, 0, 1));

    small.states[22] = NUCLEUS;
    EXPECT_TRUE(small.at(0, 0, 1) == NUCLEUS);

    small.states[26] = POTENTIAL_PORE;
    EXPECT_TRUE(small.at(26) == POTENTIAL_PORE);

    small.states = std::vector<BoxState>(30);
    small.height = 3;
    small.depth = 2;
    small.width = 5;
    small.boxes = 30;
    small.min = Vec<int>(0, 0, 0);
    small.max = Vec<int>(5, 3, 2);

    small.states[20] = POTENTIAL_PORE;
    EXPECT_TRUE(small.at(0, 1, 1) == POTENTIAL_PORE);
}

TEST_F(ProteinGridTest, add_atom) {
    Vec<int> box = Vec<int>(0, 1, 1);
    small.occupied[small.index(box)].insert(atom_1->id);
    EXPECT_EQ(1, small.occupied[small.index(box)].size());
    EXPECT_FALSE(small.at(box) == OCCUPIED);
    small.occupied[small.index(box)].insert(atom_2->id);
    EXPECT_EQ(1, small.occupied[small.index(box)].size());
    small.occupied[small.index(box)].insert(atom_3->id);
    EXPECT_EQ(2, small.occupied[small.index(box)].size());
}

TEST_F(ProteinGridTest, trace) {
    small.trace = std::vector<std::unordered_set<size_t>>(small.size());
    size_t index = small.index(0, 1, 1);
    Vec<double> vec = Vec<double>(1.0, 0.0, 0.0);
    std::shared_ptr<AtomPair> pair_1 = std::make_shared<AtomPair>(AtomPair(0, atom_1, atom_3, vec, vec.unit()));
    small.trace[index].insert(pair_1->id);
    EXPECT_EQ(1, small.trace[index].size());
    std::shared_ptr<AtomPair> pair_2 = std::make_shared<AtomPair>(AtomPair(0, atom_1, atom_3, vec, vec.unit()));
    small.trace[index].insert(pair_2->id);
    EXPECT_EQ(1, small.trace[index].size());
    std::shared_ptr<AtomPair> pair_3 = std::make_shared<AtomPair>(AtomPair(1, atom_1, atom_3, vec, vec.unit()));
    small.trace[index].insert(pair_3->id);
    EXPECT_EQ(2, small.trace[index].size());
}

TEST_F(ProteinGridTest, constructor_small) {
    // solvent radius = 1.0
    // max_vdW = 1.8
    // => 2 * margin = 2 * (1.8 + 2.0) = 7.6
    EXPECT_EQ(1.0, small.box_length);
    EXPECT_EQ(0.5, small.box_radius);
    EXPECT_EQ(1.0, small.box_volume);
    EXPECT_EQ(1.0, small.solvent_radius);
    EXPECT_EQ(3, small.atoms.size());
    EXPECT_EQ(Vec<int>(-7, -5, -6), small.min);
    EXPECT_EQ(Vec<int>(7, 8, 7), small.max);
    EXPECT_EQ(14, small.width);
    EXPECT_EQ(13, small.height);
    EXPECT_EQ(13, small.depth);
    EXPECT_EQ(2366, small.occupied.size());
    EXPECT_EQ(2366, small.states.size());
    EXPECT_EQ(0, small.trace.size());
    EXPECT_EQ(3.8, small.margin);

    EXPECT_NEAR(-3, small.atoms[0]->coord.x, 0.000001);
    EXPECT_NEAR(4.2, small.atoms[0]->coord.y, 0.000001);
    EXPECT_NEAR(1.1, small.atoms[0]->coord.z, 0.000001);

    EXPECT_NEAR(1.2, small.atoms[1]->coord.x, 0.000001);
    EXPECT_NEAR(-0.9, small.atoms[1]->coord.y, 0.000001);
    EXPECT_NEAR(3.2, small.atoms[1]->coord.z, 0.000001);

    EXPECT_NEAR(2.7, small.atoms[2]->coord.x, 0.000001);
    EXPECT_NEAR(1.3, small.atoms[2]->coord.y, 0.000001);
    EXPECT_NEAR(-2.1, small.atoms[2]->coord.z, 0.000001);

    EXPECT_EQ(C, small.atoms[0]->type);
    EXPECT_EQ(H, small.atoms[1]->type);
    EXPECT_EQ(O, small.atoms[2]->type);

    EXPECT_NEAR(0.25598000714, small.centre_of_mass.x, 0.001);
    EXPECT_NEAR(2.5434487683, small.centre_of_mass.y, 0.001);
    EXPECT_NEAR(-0.7279186005, small.centre_of_mass.z, 0.001);

    for (int x = small.min.x; x < small.max.x; x++) {
        for (int y = small.min.y; y < small.max.y; y++) {
            for (int z = small.min.z; z < small.max.z; z++) {
                Vec<int> box = Vec<int>(x, y, z);
                if (box == Vec<int>(-5, 4, 0) || box == Vec<int>(-5, 4, 1) || box == Vec<int>(-4, 3, 0)
                    || box == Vec<int>(-4, 3, 1) || box == Vec<int>(-4, 3, 2) || box == Vec<int>(-4, 4, 0)
                    || box == Vec<int>(-4, 4, 1) || box == Vec<int>(-4, 4, 2) || box == Vec<int>(-4, 5, 0)
                    || box == Vec<int>(-4, 5, 1) || box == Vec<int>(-3, 3, 0) || box == Vec<int>(-3, 3, 1)
                    || box == Vec<int>(-3, 3, 2) || box == Vec<int>(-3, 4, 0) || box == Vec<int>(-3, 4, 1)
                    || box == Vec<int>(-3, 4, 2) || box == Vec<int>(-3, 5, 0) || box == Vec<int>(-3, 5, 1)
                    || box == Vec<int>(-2, 4, 0) || box == Vec<int>(-2, 4, 1) || box == Vec<int>(0, -2, 2)
                    || box == Vec<int>(0, -2, 3) || box == Vec<int>(0, -1, 2) || box == Vec<int>(0, -1, 3)
                    || box == Vec<int>(1, -2, 2) || box == Vec<int>(1, -2, 3) || box == Vec<int>(1, -1, 2)
                    || box == Vec<int>(1, -1, 3) || box == Vec<int>(1, 0, -3) || box == Vec<int>(1, 1, -3)
                    || box == Vec<int>(1, 1, -2) || box == Vec<int>(2, 0, -3) || box == Vec<int>(2, 0, -2)
                    || box == Vec<int>(2, 1, -4) || box == Vec<int>(2, 1, -3) || box == Vec<int>(2, 1, -2)
                    || box == Vec<int>(2, 2, -3) || box == Vec<int>(2, 2, -2) || box == Vec<int>(3, 0, -3)
                    || box == Vec<int>(3, 0, -2) || box == Vec<int>(3, 1, -3) || box == Vec<int>(3, 1, -2)
                    || box == Vec<int>(3, 2, -3)) {
                    EXPECT_TRUE(small.at(box) == OCCUPIED);
                    EXPECT_EQ(1, small.occupied[small.index(box)].size());
                } else {
                    EXPECT_TRUE(small.at(box) == UNCLASSIFIED);
                    EXPECT_EQ(0, small.occupied[small.index(box)].size());
                }
            }
        }
    }
}

TEST_F(ProteinGridTest, constructor_small_05) {
    // solvent radius = 1.0
    // max_vdW = 1.8
    // => 2 * margin = 2 * (1.8 + 2.0) = 7.6
    EXPECT_EQ(0.5, small_05.box_length);
    EXPECT_EQ(0.25, small_05.box_radius);
    EXPECT_EQ(0.125, small_05.box_volume);
    EXPECT_EQ(1.0, small_05.solvent_radius);
    EXPECT_EQ(3, small_05.atoms.size());
    EXPECT_EQ(Vec<int>(-14, -10, -12), small_05.min);
    EXPECT_EQ(Vec<int>(13, 16, 14), small_05.max);
    EXPECT_EQ(26, small_05.height);
    EXPECT_EQ(27, small_05.width);
    EXPECT_EQ(26, small_05.depth);
    EXPECT_EQ(18252, small_05.occupied.size());
    EXPECT_EQ(18252, small_05.states.size());
    EXPECT_EQ(0, small_05.trace.size());

    EXPECT_NEAR(-3, small_05.atoms[0]->coord.x, 0.000001);
    EXPECT_NEAR(4.2, small_05.atoms[0]->coord.y, 0.000001);
    EXPECT_NEAR(1.1, small_05.atoms[0]->coord.z, 0.000001);

    EXPECT_NEAR(1.2, small_05.atoms[1]->coord.x, 0.000001);
    EXPECT_NEAR(-0.9, small_05.atoms[1]->coord.y, 0.000001);
    EXPECT_NEAR(3.2, small_05.atoms[1]->coord.z, 0.000001);

    EXPECT_NEAR(2.7, small_05.atoms[2]->coord.x, 0.000001);
    EXPECT_NEAR(1.3, small_05.atoms[2]->coord.y, 0.000001);
    EXPECT_NEAR(-2.1, small_05.atoms[2]->coord.z, 0.000001);

    for (int x = small_05.min.x; x < small_05.max.x; x++) {
        for (int y = small_05.min.y; y < small_05.max.y; y++) {
            for (int z = small_05.min.z; z < small_05.max.z; z++) {
                Vec<int> box = Vec<int>(x, y, z);
                if (box == Vec<int>(-9, 6, 1) || box == Vec<int>(-9, 6, 2) || box == Vec<int>(-9, 6, 3)
                    || box == Vec<int>(-9, 7, 0) || box == Vec<int>(-9, 7, 1) || box == Vec<int>(-9, 7, 2)
                    || box == Vec<int>(-9, 7, 3) || box == Vec<int>(-9, 8, 0) || box == Vec<int>(-9, 8, 1)
                    || box == Vec<int>(-9, 8, 2) || box == Vec<int>(-9, 8, 3) || box == Vec<int>(-9, 8, 4)
                    || box == Vec<int>(-9, 9, 0) || box == Vec<int>(-9, 9, 1) || box == Vec<int>(-9, 9, 2)
                    || box == Vec<int>(-9, 9, 3) || box == Vec<int>(-9, 10, 1) || box == Vec<int>(-9, 10, 2)
                    || box == Vec<int>(-8, 5, 1) || box == Vec<int>(-8, 5, 2) || box == Vec<int>(-8, 6, 0)
                    || box == Vec<int>(-8, 6, 1) || box == Vec<int>(-8, 6, 2) || box == Vec<int>(-8, 6, 3)
                    || box == Vec<int>(-8, 6, 4) || box == Vec<int>(-8, 7, -1) || box == Vec<int>(-8, 7, 0)
                    || box == Vec<int>(-8, 7, 1) || box == Vec<int>(-8, 7, 2) || box == Vec<int>(-8, 7, 3)
                    || box == Vec<int>(-8, 7, 4) || box == Vec<int>(-8, 8, -1) || box == Vec<int>(-8, 8, 0)
                    || box == Vec<int>(-8, 8, 1) || box == Vec<int>(-8, 8, 2) || box == Vec<int>(-8, 8, 3)
                    || box == Vec<int>(-8, 8, 4) || box == Vec<int>(-8, 9, -1) || box == Vec<int>(-8, 9, 0)
                    || box == Vec<int>(-8, 9, 1) || box == Vec<int>(-8, 9, 2) || box == Vec<int>(-8, 9, 3)
                    || box == Vec<int>(-8, 9, 4) || box == Vec<int>(-8, 10, 0) || box == Vec<int>(-8, 10, 1)
                    || box == Vec<int>(-8, 10, 2) || box == Vec<int>(-8, 10, 3) || box == Vec<int>(-7, 5, 0)
                    || box == Vec<int>(-7, 5, 1) || box == Vec<int>(-7, 5, 2) || box == Vec<int>(-7, 5, 3)
                    || box == Vec<int>(-7, 6, -1) || box == Vec<int>(-7, 6, 0) || box == Vec<int>(-7, 6, 1)
                    || box == Vec<int>(-7, 6, 2) || box == Vec<int>(-7, 6, 3) || box == Vec<int>(-7, 6, 4)
                    || box == Vec<int>(-7, 7, -1) || box == Vec<int>(-7, 7, 0) || box == Vec<int>(-7, 7, 1)
                    || box == Vec<int>(-7, 7, 2) || box == Vec<int>(-7, 7, 3) || box == Vec<int>(-7, 7, 4)
                    || box == Vec<int>(-7, 8, -1) || box == Vec<int>(-7, 8, 0) || box == Vec<int>(-7, 8, 1)
                    || box == Vec<int>(-7, 8, 2) || box == Vec<int>(-7, 8, 3) || box == Vec<int>(-7, 8, 4)
                    || box == Vec<int>(-7, 8, 5) || box == Vec<int>(-7, 9, -1) || box == Vec<int>(-7, 9, 0)
                    || box == Vec<int>(-7, 9, 1) || box == Vec<int>(-7, 9, 2) || box == Vec<int>(-7, 9, 3)
                    || box == Vec<int>(-7, 9, 4) || box == Vec<int>(-7, 10, 0) || box == Vec<int>(-7, 10, 1)
                    || box == Vec<int>(-7, 10, 2) || box == Vec<int>(-7, 10, 3) || box == Vec<int>(-7, 10, 4)
                    || box == Vec<int>(-7, 11, 1) || box == Vec<int>(-7, 11, 2) || box == Vec<int>(-7, 11, 3)
                    || box == Vec<int>(-6, 5, 0) || box == Vec<int>(-6, 5, 1) || box == Vec<int>(-6, 5, 2)
                    || box == Vec<int>(-6, 5, 3) || box == Vec<int>(-6, 6, -1) || box == Vec<int>(-6, 6, 0)
                    || box == Vec<int>(-6, 6, 1) || box == Vec<int>(-6, 6, 2) || box == Vec<int>(-6, 6, 3)
                    || box == Vec<int>(-6, 6, 4) || box == Vec<int>(-6, 7, -1) || box == Vec<int>(-6, 7, 0)
                    || box == Vec<int>(-6, 7, 1) || box == Vec<int>(-6, 7, 2) || box == Vec<int>(-6, 7, 3)
                    || box == Vec<int>(-6, 7, 4) || box == Vec<int>(-6, 8, -1) || box == Vec<int>(-6, 8, 0)
                    || box == Vec<int>(-6, 8, 1) || box == Vec<int>(-6, 8, 2) || box == Vec<int>(-6, 8, 3)
                    || box == Vec<int>(-6, 8, 4) || box == Vec<int>(-6, 8, 5) || box == Vec<int>(-6, 9, -1)
                    || box == Vec<int>(-6, 9, 0) || box == Vec<int>(-6, 9, 1) || box == Vec<int>(-6, 9, 2)
                    || box == Vec<int>(-6, 9, 3) || box == Vec<int>(-6, 9, 4) || box == Vec<int>(-6, 10, 0)
                    || box == Vec<int>(-6, 10, 1) || box == Vec<int>(-6, 10, 2) || box == Vec<int>(-6, 10, 3)
                    || box == Vec<int>(-6, 10, 4) || box == Vec<int>(-6, 11, 1) || box == Vec<int>(-6, 11, 2)
                    || box == Vec<int>(-6, 11, 3) || box == Vec<int>(-5, 5, 1) || box == Vec<int>(-5, 5, 2)
                    || box == Vec<int>(-5, 6, 0) || box == Vec<int>(-5, 6, 1) || box == Vec<int>(-5, 6, 2)
                    || box == Vec<int>(-5, 6, 3) || box == Vec<int>(-5, 6, 4) || box == Vec<int>(-5, 7, -1)
                    || box == Vec<int>(-5, 7, 0) || box == Vec<int>(-5, 7, 1) || box == Vec<int>(-5, 7, 2)
                    || box == Vec<int>(-5, 7, 3) || box == Vec<int>(-5, 7, 4) || box == Vec<int>(-5, 8, -1)
                    || box == Vec<int>(-5, 8, 0) || box == Vec<int>(-5, 8, 1) || box == Vec<int>(-5, 8, 2)
                    || box == Vec<int>(-5, 8, 3) || box == Vec<int>(-5, 8, 4) || box == Vec<int>(-5, 9, -1)
                    || box == Vec<int>(-5, 9, 0) || box == Vec<int>(-5, 9, 1) || box == Vec<int>(-5, 9, 2)
                    || box == Vec<int>(-5, 9, 3) || box == Vec<int>(-5, 9, 4) || box == Vec<int>(-5, 10, 0)
                    || box == Vec<int>(-5, 10, 1) || box == Vec<int>(-5, 10, 2) || box == Vec<int>(-5, 10, 3)
                    || box == Vec<int>(-4, 6, 1) || box == Vec<int>(-4, 6, 2) || box == Vec<int>(-4, 6, 3)
                    || box == Vec<int>(-4, 7, 0) || box == Vec<int>(-4, 7, 1) || box == Vec<int>(-4, 7, 2)
                    || box == Vec<int>(-4, 7, 3) || box == Vec<int>(-4, 8, 0) || box == Vec<int>(-4, 8, 1)
                    || box == Vec<int>(-4, 8, 2) || box == Vec<int>(-4, 8, 3) || box == Vec<int>(-4, 8, 4)
                    || box == Vec<int>(-4, 9, 0) || box == Vec<int>(-4, 9, 1) || box == Vec<int>(-4, 9, 2)
                    || box == Vec<int>(-4, 9, 3) || box == Vec<int>(-4, 10, 1) || box == Vec<int>(-4, 10, 2)
                    || box == Vec<int>(0, -3, 5) || box == Vec<int>(0, -3, 6) || box == Vec<int>(0, -3, 7)
                    || box == Vec<int>(0, -2, 5) || box == Vec<int>(0, -2, 6) || box == Vec<int>(0, -2, 7)
                    || box == Vec<int>(0, -1, 6) || box == Vec<int>(1, -4, 5) || box == Vec<int>(1, -4, 6)
                    || box == Vec<int>(1, -4, 7) || box == Vec<int>(1, -3, 4) || box == Vec<int>(1, -3, 5)
                    || box == Vec<int>(1, -3, 6) || box == Vec<int>(1, -3, 7) || box == Vec<int>(1, -3, 8)
                    || box == Vec<int>(1, -2, 4) || box == Vec<int>(1, -2, 5) || box == Vec<int>(1, -2, 6)
                    || box == Vec<int>(1, -2, 7) || box == Vec<int>(1, -2, 8) || box == Vec<int>(1, -1, 5)
                    || box == Vec<int>(1, -1, 6) || box == Vec<int>(1, -1, 7) || box == Vec<int>(2, -4, 5)
                    || box == Vec<int>(2, -4, 6) || box == Vec<int>(2, -4, 7) || box == Vec<int>(2, -3, 4)
                    || box == Vec<int>(2, -3, 5) || box == Vec<int>(2, -3, 6) || box == Vec<int>(2, -3, 7)
                    || box == Vec<int>(2, -3, 8) || box == Vec<int>(2, -2, 4) || box == Vec<int>(2, -2, 5)
                    || box == Vec<int>(2, -2, 6) || box == Vec<int>(2, -2, 7) || box == Vec<int>(2, -2, 8)
                    || box == Vec<int>(2, -1, 4) || box == Vec<int>(2, -1, 5) || box == Vec<int>(2, -1, 6)
                    || box == Vec<int>(2, -1, 7) || box == Vec<int>(2, 0, 6) || box == Vec<int>(2, 2, -5)
                    || box == Vec<int>(2, 2, -4) || box == Vec<int>(3, -4, 5) || box == Vec<int>(3, -4, 6)
                    || box == Vec<int>(3, -4, 7) || box == Vec<int>(3, -3, 4) || box == Vec<int>(3, -3, 5)
                    || box == Vec<int>(3, -3, 6) || box == Vec<int>(3, -3, 7) || box == Vec<int>(3, -2, 4)
                    || box == Vec<int>(3, -2, 5) || box == Vec<int>(3, -2, 6) || box == Vec<int>(3, -2, 7)
                    || box == Vec<int>(3, -2, 8) || box == Vec<int>(3, -1, 5) || box == Vec<int>(3, -1, 6)
                    || box == Vec<int>(3, -1, 7) || box == Vec<int>(3, 0, -5) || box == Vec<int>(3, 0, -4)
                    || box == Vec<int>(3, 1, -6) || box == Vec<int>(3, 1, -5) || box == Vec<int>(3, 1, -4)
                    || box == Vec<int>(3, 1, -3) || box == Vec<int>(3, 2, -7) || box == Vec<int>(3, 2, -6)
                    || box == Vec<int>(3, 2, -5) || box == Vec<int>(3, 2, -4) || box == Vec<int>(3, 2, -3)
                    || box == Vec<int>(3, 3, -6) || box == Vec<int>(3, 3, -5) || box == Vec<int>(3, 3, -4)
                    || box == Vec<int>(3, 3, -3) || box == Vec<int>(3, 4, -6) || box == Vec<int>(3, 4, -5)
                    || box == Vec<int>(3, 4, -4) || box == Vec<int>(4, -3, 5) || box == Vec<int>(4, -3, 6)
                    || box == Vec<int>(4, -2, 5) || box == Vec<int>(4, -2, 6) || box == Vec<int>(4, -2, 7)
                    || box == Vec<int>(4, 0, -6) || box == Vec<int>(4, 0, -5) || box == Vec<int>(4, 0, -4)
                    || box == Vec<int>(4, 0, -3) || box == Vec<int>(4, 1, -7) || box == Vec<int>(4, 1, -6)
                    || box == Vec<int>(4, 1, -5) || box == Vec<int>(4, 1, -4) || box == Vec<int>(4, 1, -3)
                    || box == Vec<int>(4, 2, -7) || box == Vec<int>(4, 2, -6) || box == Vec<int>(4, 2, -5)
                    || box == Vec<int>(4, 2, -4) || box == Vec<int>(4, 2, -3) || box == Vec<int>(4, 2, -2)
                    || box == Vec<int>(4, 3, -7) || box == Vec<int>(4, 3, -6) || box == Vec<int>(4, 3, -5)
                    || box == Vec<int>(4, 3, -4) || box == Vec<int>(4, 3, -3) || box == Vec<int>(4, 3, -2)
                    || box == Vec<int>(4, 4, -6) || box == Vec<int>(4, 4, -5) || box == Vec<int>(4, 4, -4)
                    || box == Vec<int>(4, 4, -3) || box == Vec<int>(5, 0, -6) || box == Vec<int>(5, 0, -5)
                    || box == Vec<int>(5, 0, -4) || box == Vec<int>(5, 0, -3) || box == Vec<int>(5, 1, -7)
                    || box == Vec<int>(5, 1, -6) || box == Vec<int>(5, 1, -5) || box == Vec<int>(5, 1, -4)
                    || box == Vec<int>(5, 1, -3) || box == Vec<int>(5, 1, -2) || box == Vec<int>(5, 2, -7)
                    || box == Vec<int>(5, 2, -6) || box == Vec<int>(5, 2, -5) || box == Vec<int>(5, 2, -4)
                    || box == Vec<int>(5, 2, -3) || box == Vec<int>(5, 2, -2) || box == Vec<int>(5, 3, -7)
                    || box == Vec<int>(5, 3, -6) || box == Vec<int>(5, 3, -5) || box == Vec<int>(5, 3, -4)
                    || box == Vec<int>(5, 3, -3) || box == Vec<int>(5, 3, -2) || box == Vec<int>(5, 4, -7)
                    || box == Vec<int>(5, 4, -6) || box == Vec<int>(5, 4, -5) || box == Vec<int>(5, 4, -4)
                    || box == Vec<int>(5, 4, -3) || box == Vec<int>(5, 5, -5) || box == Vec<int>(5, 5, -4)
                    || box == Vec<int>(6, 0, -6) || box == Vec<int>(6, 0, -5) || box == Vec<int>(6, 0, -4)
                    || box == Vec<int>(6, 0, -3) || box == Vec<int>(6, 1, -7) || box == Vec<int>(6, 1, -6)
                    || box == Vec<int>(6, 1, -5) || box == Vec<int>(6, 1, -4) || box == Vec<int>(6, 1, -3)
                    || box == Vec<int>(6, 2, -7) || box == Vec<int>(6, 2, -6) || box == Vec<int>(6, 2, -5)
                    || box == Vec<int>(6, 2, -4) || box == Vec<int>(6, 2, -3) || box == Vec<int>(6, 2, -2)
                    || box == Vec<int>(6, 3, -7) || box == Vec<int>(6, 3, -6) || box == Vec<int>(6, 3, -5)
                    || box == Vec<int>(6, 3, -4) || box == Vec<int>(6, 3, -3) || box == Vec<int>(6, 4, -6)
                    || box == Vec<int>(6, 4, -5) || box == Vec<int>(6, 4, -4) || box == Vec<int>(6, 4, -3)
                    || box == Vec<int>(7, 0, -5) || box == Vec<int>(7, 1, -6) || box == Vec<int>(7, 1, -5)
                    || box == Vec<int>(7, 1, -4) || box == Vec<int>(7, 1, -3) || box == Vec<int>(7, 2, -6)
                    || box == Vec<int>(7, 2, -5) || box == Vec<int>(7, 2, -4) || box == Vec<int>(7, 2, -3)
                    || box == Vec<int>(7, 3, -6) || box == Vec<int>(7, 3, -5) || box == Vec<int>(7, 3, -4)
                    || box == Vec<int>(7, 3, -3) || box == Vec<int>(7, 4, -5) || box == Vec<int>(7, 4, -4)) {
                    EXPECT_TRUE(small_05.at(box) == OCCUPIED);
                    EXPECT_EQ(1, small_05.occupied[small_05.index(box)].size());
                } else {
                    EXPECT_TRUE(small_05.at(box) == UNCLASSIFIED);
                    EXPECT_EQ(0, small_05.occupied[small_05.index(box)].size());
                }
            }
        }
    }
}

TEST_F(ProteinGridTest, constructor_g_per) {
    // solvent radius = 1.0
    // max_vdW = 1.8
    // => 2 * margin = 2 * (1.8 + 2.0) = 7.6
    EXPECT_EQ(1.0, g_per.box_length);
    EXPECT_EQ(0.5, g_per.box_radius);
    EXPECT_EQ(1, g_per.box_volume);
    EXPECT_EQ(1.0, g_per.solvent_radius);
    EXPECT_EQ(7, g_per.atoms.size());
    EXPECT_EQ(Vec<int>(0, 0, 0), g_per.min);
    EXPECT_EQ(Vec<int>(24, 38, 8), g_per.max);
    EXPECT_EQ(38, g_per.height);
    EXPECT_EQ(24, g_per.width);
    EXPECT_EQ(8, g_per.depth);
    EXPECT_EQ(7296, g_per.occupied.size());
    EXPECT_EQ(7296, g_per.states.size());
    EXPECT_EQ(0, g_per.trace.size());

    for (int x = g_per.min.x; x < g_per.max.x; x++) {
        for (int y = g_per.min.y; y < g_per.max.y; y++) {
            for (int z = g_per.min.z; z < g_per.max.z; z++) {
                Vec<int> box = Vec<int>(x, y, z);
                if (box == Vec<int>(3, 18, 3) || box == Vec<int>(3, 18, 4) || box == Vec<int>(3, 19, 3)
                    || box == Vec<int>(3, 19, 4) || box == Vec<int>(4, 18, 3) || box == Vec<int>(4, 18, 4)
                    || box == Vec<int>(4, 19, 3) || box == Vec<int>(6, 33, 3) || box == Vec<int>(6, 33, 4)
                    || box == Vec<int>(6, 34, 3) || box == Vec<int>(7, 3, 3) || box == Vec<int>(7, 3, 4)
                    || box == Vec<int>(7, 4, 3) || box == Vec<int>(7, 33, 3) || box == Vec<int>(7, 33, 4)
                    || box == Vec<int>(7, 34, 3) || box == Vec<int>(7, 34, 4) || box == Vec<int>(8, 3, 3)
                    || box == Vec<int>(8, 3, 4) || box == Vec<int>(8, 4, 3) || box == Vec<int>(8, 4, 4)
                    || box == Vec<int>(9, 3, 3) || box == Vec<int>(11, 18, 3) || box == Vec<int>(11, 18, 4)
                    || box == Vec<int>(11, 19, 3) || box == Vec<int>(11, 19, 4) || box == Vec<int>(12, 18, 3)
                    || box == Vec<int>(12, 18, 4) || box == Vec<int>(12, 19, 3) || box == Vec<int>(13, 33, 3)
                    || box == Vec<int>(13, 33, 4) || box == Vec<int>(13, 34, 3) || box == Vec<int>(13, 34, 4)
                    || box == Vec<int>(14, 33, 3) || box == Vec<int>(14, 33, 4) || box == Vec<int>(14, 34, 3)
                    || box == Vec<int>(17, 3, 3) || box == Vec<int>(17, 3, 4) || box == Vec<int>(17, 4, 3)
                    || box == Vec<int>(17, 4, 4) || box == Vec<int>(18, 3, 3) || box == Vec<int>(18, 3, 4)
                    || box == Vec<int>(18, 4, 3) || box == Vec<int>(19, 18, 3) || box == Vec<int>(19, 18, 4)
                    || box == Vec<int>(19, 19, 3) || box == Vec<int>(19, 19, 4) || box == Vec<int>(20, 18, 3)
                    || box == Vec<int>(20, 18, 4) || box == Vec<int>(20, 19, 3)) {
                    EXPECT_TRUE(g_per.at(box) == OCCUPIED);
                    EXPECT_EQ(1, g_per.occupied[g_per.index(box)].size());
                } else {
                    EXPECT_TRUE(g_per.at(box) == UNCLASSIFIED);
                    EXPECT_EQ(0, g_per.occupied[g_per.index(box)].size());
                }
            }
        }
    }
}

TEST_F(ProteinGridTest, build_kd_tree) {
    EXPECT_FALSE(kd.KD_constructed);
    EXPECT_EQ(nullptr, kd.KD_root);
    EXPECT_EQ(15, kd.atoms.size());

    kd.build_KD_tree();

    EXPECT_EQ(4, kd.KD_root->id);
    EXPECT_TRUE(kd.KD_root->has_left());
    EXPECT_EQ(5, kd.KD_root->KD_left->id);
    EXPECT_EQ(11, kd.KD_root->KD_left->KD_left->id);
    EXPECT_EQ(3, kd.KD_root->KD_left->KD_left->KD_left->id);
    EXPECT_EQ(14, kd.KD_root->KD_left->KD_left->KD_right->id);
    EXPECT_EQ(9, kd.KD_root->KD_left->KD_right->id);
    EXPECT_EQ(2, kd.KD_root->KD_left->KD_right->KD_left->id);
    EXPECT_EQ(8, kd.KD_root->KD_left->KD_right->KD_right->id);
    EXPECT_EQ(0, kd.KD_root->KD_right->id);
    EXPECT_EQ(10, kd.KD_root->KD_right->KD_left->id);
    EXPECT_EQ(7, kd.KD_root->KD_right->KD_left->KD_left->id);
    EXPECT_EQ(12, kd.KD_root->KD_right->KD_left->KD_right->id);
    EXPECT_EQ(13, kd.KD_root->KD_right->KD_right->id);
    EXPECT_EQ(1, kd.KD_root->KD_right->KD_right->KD_left->id);
    EXPECT_EQ(6, kd.KD_root->KD_right->KD_right->KD_right->id);
}

TEST_F(ProteinGridTest, KD_closest) {
    kd.build_KD_tree();
    Vec<double> query = Vec<double>(36.5, 10.5, 48.5);
    EXPECT_NEAR(1.0516, std::get<0>(kd.closest(query)), 0.001);

    query = Vec<double>(37.8, 8.8, 47.8);
    EXPECT_NEAR(-1.2, std::get<0>(kd.closest(query)), 0.01);

    query = Vec<double>(100, 100, 100);
    EXPECT_NEAR(90.61, std::get<0>(kd.closest(query)), 0.01);

    query = Vec<double>(3.8, 3.8, 3.8);
    EXPECT_NEAR(36.363, std::get<0>(kd.closest(query)), 0.01);

    EXPECT_NEAR(36.882, std::get<0>(kd.closest(3, 3, 3)), 0.01);

    query = Vec<double>(37.801, 22.8, 49.8);
    EXPECT_NEAR(2.801, std::get<0>(kd.closest(query)), 0.01);
}

TEST_F(ProteinGridTest, add_atom_pair) {
    small.add_atom_pair(atom_1, atom_2, Vec<double>(1.0, 2.0, -1.0), Vec<double>(-1.0, 0.0, 0.0));
    small.add_atom_pair(atom_2, atom_3, Vec<double>(3.0, 4.0, 1.0), Vec<double>(0.0, -1.0, 0.0));
    EXPECT_EQ(2, small.pairs.size());
    EXPECT_EQ(0, small.pairs[0]->id);
    EXPECT_EQ(0, small.pairs[0]->atom_1->id);
    EXPECT_EQ(0, small.pairs[0]->atom_2->id);
    EXPECT_EQ(Vec<double>(1.0, 2.0, -1.0), small.pairs[0]->vec);
    EXPECT_EQ(Vec<double>(-1.0, 0.0, 0.0), small.pairs[0]->unit);
    EXPECT_EQ(1, small.pairs[1]->id);
    EXPECT_EQ(0, small.pairs[1]->atom_1->id);
    EXPECT_EQ(1, small.pairs[1]->atom_2->id);
    EXPECT_EQ(Vec<double>(3.0, 4.0, 1.0), small.pairs[1]->vec);
    EXPECT_EQ(Vec<double>(0.0, -1.0, 0.0), small.pairs[1]->unit);
}

TEST_F(ProteinGridTest, new_cluster) {
    small.new_cluster();
    small.new_cluster();
    small.new_cluster();
    EXPECT_EQ(3, small.clusters.size());
    EXPECT_EQ(0, small.clusters[0].id);
    EXPECT_EQ(1, small.clusters[1].id);
    EXPECT_EQ(2, small.clusters[2].id);
    small.clusters[0].pore = true;
    small.clusters[2].pore = true;
    small.clusters[1].pore = false;
    EXPECT_EQ("0_pore", small.clusters[0].name());
    EXPECT_EQ("2_pore", small.clusters[2].name());
    EXPECT_EQ("1_cavity", small.clusters[1].name());

    for (size_t i = 0; i < 30; i++) {
        small.clusters[1].add_box(i, Vec<int>(0, 0, 0));
        small.clusters[2].add_box(i, Vec<int>(0, 0, 0));
    }

    EXPECT_EQ(30, small.clusters[1].size());
    EXPECT_EQ(30, small.clusters[2].size());
    EXPECT_EQ(0, small.clusters[0].size());
    for (size_t i = 0; i < 30; i++) {
        EXPECT_TRUE(small.clusters[1].has_index(i));
        EXPECT_TRUE(small.clusters[2].has_index(i));
        EXPECT_FALSE(small.clusters[0].has_index(i));
        small.clusters[1].add_box(i, Vec<int>(0, 0, 0));
        small.clusters[2].add_box(i, Vec<int>(0, 0, 0));
    }
    EXPECT_EQ(30, small.clusters[1].size());
    EXPECT_EQ(30, small.clusters[2].size());
    EXPECT_EQ(0, small.clusters[0].size());

    small.filter_pores(30, REMOVE_NOTHING);
    EXPECT_EQ(2, small.clusters.size());
    EXPECT_EQ(0, small.clusters[0].id);
    EXPECT_EQ(1, small.clusters[1].id);
    EXPECT_FALSE(small.clusters[0].pore);
    EXPECT_TRUE(small.clusters[1].pore);
    small.filter_pores(30, REMOVE_CAVITIES);
    EXPECT_EQ(1, small.clusters.size());
    EXPECT_EQ(0, small.clusters[0].id);
    small.filter_pores(30, REMOVE_PORES);
    EXPECT_EQ(0, small.clusters.size());
    EXPECT_EQ(0, small.next_cluster_id);
}