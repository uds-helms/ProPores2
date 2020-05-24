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
#include <vector>
#include "atom.h"
#include "grid.h"
#include "enums.h"
#include "reader.h"
#include "vector.h"
#include "pore_id.h"
#include "gtest/gtest.h"


class PoreIDTest : public ::testing::Test {
public:
    PoreIDTest() :
            g_psp("test_files/psp_scan.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true),
            a_psp(g_psp.atoms[0]),
            g_col("test_files/collision_detection.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true),
            g_cyl("test_files/cylinder_completion.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true),
            g_per("test_files/perpendicular.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true) {
        g_psp.states = std::vector<BoxState>(g_psp.size());

        // perpendicular
        psp_scan(g_per);
        collision_detection(g_per);
        std::vector<std::shared_ptr<AtomPair>> pairs;
        size_t id = 0;
        for (const std::shared_ptr<AtomPair> &p: g_per.pairs) {
            if ((p->atom_1->id == 0 && p->atom_2->id == 1) || (p->atom_1->id == 1 && p->atom_2->id == 2)
                || (p->atom_1->id == 3 && p->atom_2->id == 4) || (p->atom_1->id == 5 && p->atom_2->id == 6)) {
                pairs.emplace_back(new AtomPair(id, p->atom_1, p->atom_2, p->vec, p->unit));
                id++;
            }
        }

        g_per.pairs.clear();
        for (const std::shared_ptr<AtomPair> &p: pairs) {
            g_per.pairs.emplace_back(new AtomPair(p->id, p->atom_1, p->atom_2, p->vec, p->unit));
        }

        cylinder_completion(g_per);
        perpendicularity_cylinder(g_per);
    }

    // grid of 8x8x8 for PSP scan tests
    ProteinGrid g_psp;
    std::shared_ptr<Atom> &a_psp;
    // grid for collision detection tests
    ProteinGrid g_col;
    // grid for cylinder completion tests
    ProteinGrid g_cyl;
    // grid for perpendicularity tests
    ProteinGrid g_per;
};


TEST_F(PoreIDTest, check_dimensions) {
    EXPECT_EQ(Vec<int>(8, 8, 8), g_psp.max);
    EXPECT_EQ(Vec<int>(0, 0, 0), g_psp.min);
}

TEST_F(PoreIDTest, PSP_scan_outside) {
    g_psp.at(0, 0, 0) = OCCUPIED;
    g_psp.at(7, 0, 0) = OCCUPIED;
    g_psp.at(0, 7, 0) = OCCUPIED;
    g_psp.at(0, 7, 7) = OCCUPIED;

    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (int x = g_psp.min.x; x < g_psp.max.x; x++) {
        for (int y = g_psp.min.y; y < g_psp.max.y; y++) {
            for (int z = g_psp.min.z; z < g_psp.max.z; z++) {
                EXPECT_FALSE(g_psp.at(x, y, z) == POTENTIAL_PORE);

                if (y == 0 && z == 0 && x > 0 && x < 7) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 0 && z == 0 && y > 0 && y < 7) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 0 && y == 7 && z > 0 && z < 7) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else {
                    EXPECT_EQ(0, psp[g_psp.index(x, y, z)]);
                }
            }
        }
    }
}

TEST_F(PoreIDTest, PSP_scan_outside_middle) {
    g_psp.at(0, 0, 0) = OCCUPIED;
    g_psp.at(4, 0, 0) = OCCUPIED;
    g_psp.at(7, 0, 0) = OCCUPIED;
    g_psp.at(0, 4, 0) = OCCUPIED;
    g_psp.at(0, 7, 0) = OCCUPIED;
    g_psp.at(0, 7, 7) = OCCUPIED;
    g_psp.at(0, 7, 4) = OCCUPIED;

    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (int x = g_psp.min.x; x < g_psp.max.x; x++) {
        for (int y = g_psp.min.y; y < g_psp.max.y; y++) {
            for (int z = g_psp.min.z; z < g_psp.max.z; z++) {
                EXPECT_FALSE(g_psp.at(x, y, z) == POTENTIAL_PORE);

                if (y == 0 && z == 0 && x > 0 && x < 7 && x != 4) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 0 && z == 0 && y > 0 && y < 7 && y != 4) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 0 && y == 7 && z > 0 && z < 7 && z != 4) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else {
                    EXPECT_EQ(0, psp[g_psp.index(x, y, z)]);
                }
            }
        }
    }
}

TEST_F(PoreIDTest, PSP_scan_near_outside) {
    g_psp.at(1, 1, 1) = OCCUPIED;
    g_psp.at(6, 1, 1) = OCCUPIED;
    g_psp.at(1, 6, 1) = OCCUPIED;
    g_psp.at(1, 6, 6) = OCCUPIED;

    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (int x = g_psp.min.x; x < g_psp.max.x; x++) {
        for (int y = g_psp.min.y; y < g_psp.max.y; y++) {
            for (int z = g_psp.min.z; z < g_psp.max.z; z++) {
                EXPECT_FALSE(g_psp.at(x, y, z) == POTENTIAL_PORE);
                if (y == 1 && z == 1 && x > 1 && x < 6) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 1 && z == 1 && y > 1 && y < 6) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 1 && y == 6 && z > 1 && z < 6) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else {
                    EXPECT_EQ(0, psp[g_psp.index(x, y, z)]);
                }
            }
        }
    }
}

TEST_F(PoreIDTest, PSP_scan_near_outside_2) {
    g_psp.at(1, 1, 1) = OCCUPIED;
    g_psp.at(6, 1, 1) = OCCUPIED;
    g_psp.at(1, 6, 1) = OCCUPIED;
    g_psp.at(1, 6, 6) = OCCUPIED;
    g_psp.at(7, 1, 1) = OCCUPIED;
    g_psp.at(1, 7, 1) = OCCUPIED;
    g_psp.at(1, 6, 7) = OCCUPIED;

    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (int x = g_psp.min.x; x < g_psp.max.x; x++) {
        for (int y = g_psp.min.y; y < g_psp.max.y; y++) {
            for (int z = g_psp.min.z; z < g_psp.max.z; z++) {
                EXPECT_FALSE(g_psp.at(x, y, z) == POTENTIAL_PORE);
                if (y == 1 && z == 1 && x > 1 && x < 6) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 1 && z == 1 && y > 1 && y < 6) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else if (x == 1 && y == 6 && z > 1 && z < 6) {
                    EXPECT_EQ(1, psp[g_psp.index(x, y, z)]);
                } else {
                    EXPECT_EQ(0, psp[g_psp.index(x, y, z)]);
                }
            }
        }
    }
}

TEST_F(PoreIDTest, PSP_scan_only_one_in_each) {
    g_psp.at(0, 0, 0) = OCCUPIED;
    g_psp.at(1, 1, 1) = OCCUPIED;
    g_psp.at(2, 2, 2) = OCCUPIED;

    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (size_t i = 0; i < g_psp.size(); i++) {
        EXPECT_TRUE(g_psp.states[i] != POTENTIAL_PORE);
        EXPECT_EQ(0, psp[i]);
    }
}

TEST_F(PoreIDTest, PSP_scan_empty) {
    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (size_t i = 0; i < g_psp.size(); i++) {
        EXPECT_TRUE(g_psp.states[i] != POTENTIAL_PORE);
        EXPECT_EQ(0, psp[i]);
    }
}

TEST_F(PoreIDTest, PSP_scan_one_full) {
    for (int x = g_psp.min.x; x < g_psp.max.x; x++) {
        g_psp.at(x, 1, 1) = OCCUPIED;
    }
    std::vector<uint8_t> psp = psp_scan(g_psp);
    for (size_t i = 0; i < g_psp.size(); i++) {
        EXPECT_TRUE(g_psp.states[i] != POTENTIAL_PORE);
        EXPECT_EQ(0, psp[i]);
    }
}

TEST_F(PoreIDTest, PSP_scan_multiple) {
    // surround (3, 3, 3)
    g_psp.at(2, 3, 3) = OCCUPIED;
    g_psp.at(3, 2, 3) = OCCUPIED;
    g_psp.at(3, 3, 2) = OCCUPIED;
    g_psp.at(4, 3, 3) = OCCUPIED;
    g_psp.at(3, 4, 3) = OCCUPIED;
    g_psp.at(3, 3, 4) = OCCUPIED;

    std::vector<uint8_t> psp = psp_scan(g_psp);

    for (int x = g_psp.min.x; x < g_psp.max.x; x++) {
        for (int y = g_psp.min.y; y < g_psp.max.y; y++) {
            for (int z = g_psp.min.z; z < g_psp.max.z; z++) {
                if (y == 3 && z == 3 && x == 3) {
                    EXPECT_TRUE(g_psp.at(x, y, z) == POTENTIAL_PORE);
                    EXPECT_EQ(3, psp[g_psp.index(x, y, z)]);
                } else {
                    EXPECT_TRUE(g_psp.at(x, y, z) != POTENTIAL_PORE);
                    EXPECT_EQ(0, psp[g_psp.index(x, y, z)]);
                }
            }
        }
    }
}

TEST_F(PoreIDTest, collision_empty) {
    collision_detection(g_psp);
    EXPECT_TRUE(g_psp.pairs.empty());
}

TEST_F(PoreIDTest, collision_on_grid) {
    EXPECT_EQ(5, g_col.atoms.size());
    collision_detection(g_col);

    EXPECT_EQ(7, g_col.pairs.size());
    EXPECT_EQ(0, g_col.pairs[0]->atom_1->id);
    EXPECT_EQ(1, g_col.pairs[0]->atom_2->id);
    EXPECT_EQ(0, g_col.pairs[1]->atom_1->id);
    EXPECT_EQ(4, g_col.pairs[1]->atom_2->id);
    EXPECT_EQ(1, g_col.pairs[2]->atom_1->id);
    EXPECT_EQ(2, g_col.pairs[2]->atom_2->id);
    EXPECT_EQ(1, g_col.pairs[3]->atom_1->id);
    EXPECT_EQ(3, g_col.pairs[3]->atom_2->id);
    EXPECT_EQ(1, g_col.pairs[4]->atom_1->id);
    EXPECT_EQ(4, g_col.pairs[4]->atom_2->id);
    EXPECT_EQ(2, g_col.pairs[5]->atom_1->id);
    EXPECT_EQ(3, g_col.pairs[5]->atom_2->id);
    EXPECT_EQ(3, g_col.pairs[6]->atom_1->id);
    EXPECT_EQ(4, g_col.pairs[6]->atom_2->id);
}


TEST_F(PoreIDTest, cylinder_trace) {
    psp_scan(g_cyl);
    collision_detection(g_cyl);

    std::shared_ptr<AtomPair> p = g_cyl.pairs[0];
    double r_sq = pow(1.2, 2);
    int r_grid = 2;
    double threshold = pow(r_grid - g_cyl.box_radius * sqrt(3), 2) / g_cyl.box_length_sq;
    g_cyl.trace = std::vector<std::unordered_set<size_t>>(g_cyl.size());
    cylinder_trace(g_cyl, p, r_sq, r_grid, threshold, 2.0);

    for (int x = g_cyl.min.x; x < g_cyl.max.x; x++) {
        for (int y = g_cyl.min.y; y < g_cyl.max.y; y++) {
            for (int z = g_cyl.min.z; z < g_cyl.max.z; z++) {
                Vec<int> v = Vec<int>(x, y, z);
                // (4, 3, 3) overlaps with atom 1 vdW-radius
                if (v == Vec<int>(5, 3, 3) || v == Vec<int>(6, 3, 3) || v == Vec<int>(5, 4, 3)
                    || v == Vec<int>(5, 3, 4) || v == Vec<int>(5, 4, 4) || v == Vec<int>(5, 2, 3)
                    || v == Vec<int>(5, 3, 2) || v == Vec<int>(6, 3, 4) || v == Vec<int>(6, 4, 3)) {
                    EXPECT_EQ(1, g_cyl.trace[g_cyl.index(v)].size());
                } else {
                    EXPECT_EQ(0, g_cyl.trace[g_cyl.index(v)].size());
                }
            }
        }
    }
}

TEST_F(PoreIDTest, cylinder_completion) {
    psp_scan(g_cyl);
    g_cyl.trace = std::vector<std::unordered_set<size_t>>(g_cyl.size());
    collision_detection(g_cyl);
    cylinder_completion(g_cyl);

    for (int x = g_cyl.min.x; x < g_cyl.max.x; x++) {
        for (int y = g_cyl.min.y; y < g_cyl.max.y; y++) {
            for (int z = g_cyl.min.z; z < g_cyl.max.z; z++) {
                Vec<int> v = Vec<int>(x, y, z);
                if ((4 <= x && x <= 10) || (12 <= x && x <= 18)) {
                    if ((distance(v, Vec<int>(x, 3, 3)) <= 1
                         || distance(g_cyl.centre(v), Vec<double>(x + 0.8, 3.8, 3.8)) <= 1.2)
                        && g_cyl.at(v) != OCCUPIED) {
                        EXPECT_EQ(1, g_cyl.trace[g_cyl.index(v)].size());
                    } else {
                        EXPECT_EQ(0, g_cyl.trace[g_cyl.index(v)].size());
                    }
                } else {
                    EXPECT_EQ(0, g_cyl.trace[g_cyl.index(v)].size());
                }
            }
        }
    }
}


TEST_F(PoreIDTest, perpendicular) {
    for (int x = g_per.min.x; x < g_per.max.x; x++) {
        for (int y = g_per.min.y; y < g_per.max.y; y++) {
            for (int z = g_per.min.z; z < g_per.max.z; z++) {
                Vec<int> v = Vec<int>(x, y, z);
                if (v == Vec<int>(6, 17, 3) || v == Vec<int>(6, 18, 3) || v == Vec<int>(6, 19, 3)
                    || v == Vec<int>(7, 17, 3) || v == Vec<int>(7, 18, 2) || v == Vec<int>(7, 18, 3)
                    || v == Vec<int>(7, 18, 4) || v == Vec<int>(7, 19, 3) || v == Vec<int>(7, 19, 4)
                    || v == Vec<int>(8, 17, 3) || v == Vec<int>(8, 18, 3) || v == Vec<int>(8, 18, 4)
                    || v == Vec<int>(8, 19, 3) || v == Vec<int>(8, 19, 4)) {
                    EXPECT_TRUE(g_per.at(x, y, z) == POTENTIAL_PORE);
                } else {
                    EXPECT_FALSE(g_per.at(x, y, z) == POTENTIAL_PORE);
                }
            }
        }
    }
}


TEST_F(PoreIDTest, seal_the_gaps_no_cavity) {
    seal_gaps(g_per);

    for (int x = g_per.min.x; x < g_per.max.x; x++) {
        for (int y = g_per.min.y; y < g_per.max.y; y++) {
            for (int z = g_per.min.z; z < g_per.max.z; z++) {
                Vec<int> v = Vec<int>(x, y, z);
                if (v == Vec<int>(6, 17, 3) || v == Vec<int>(6, 18, 3) || v == Vec<int>(6, 19, 3)
                    || v == Vec<int>(7, 17, 3) || v == Vec<int>(7, 18, 2) || v == Vec<int>(7, 18, 3)
                    || v == Vec<int>(7, 18, 4) || v == Vec<int>(7, 19, 3) || v == Vec<int>(7, 19, 4)
                    || v == Vec<int>(8, 17, 3) || v == Vec<int>(8, 18, 3) || v == Vec<int>(8, 18, 4)
                    || v == Vec<int>(8, 19, 3) || v == Vec<int>(8, 19, 4)) {
                    EXPECT_TRUE(g_per.at(x, y, z) == POTENTIAL_PORE);
                } else if (g_per.at(v) != OCCUPIED) {
                    EXPECT_TRUE(g_per.at(x, y, z) == BACKGROUND);
                }
            }
        }
    }
}


TEST(pore_id_tests, seal_gaps_pore) {
    ProteinGrid g_pore = ProteinGrid("test_files/pore.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_pore);
    cylinder_completion(g_pore);
    perpendicularity_cylinder(g_pore);
    std::vector<BoxState> states = std::vector<BoxState>(g_pore.states);
    seal_gaps(g_pore);

    for (int x = g_pore.min.x; x < g_pore.max.x; x++) {
        for (int y = g_pore.min.y; y < g_pore.max.y; y++) {
            for (int z = g_pore.min.z; z < g_pore.max.z; z++) {
                uint8_t state = g_pore.at(x, y, z);
                Vec<int> b = Vec<int>(x, y, z);

                if (state == OCCUPIED) {
                } else if (3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12) {
                    EXPECT_TRUE(state == POTENTIAL_PORE);

                } else if (b == Vec<int>(2, 6, 8) || b == Vec<int>(2, 7, 6) || b == Vec<int>(2, 7, 7)
                           || b == Vec<int>(2, 7, 8) || b == Vec<int>(2, 7, 9) || b == Vec<int>(2, 8, 6)
                           || b == Vec<int>(2, 8, 7) || b == Vec<int>(2, 8, 8) || b == Vec<int>(2, 8, 9)
                           || b == Vec<int>(2, 8, 10) || b == Vec<int>(2, 9, 7) || b == Vec<int>(2, 9, 8)
                           || b == Vec<int>(2, 9, 9) || b == Vec<int>(2, 10, 8) || b == Vec<int>(2, 10, 9)) {
                    EXPECT_FALSE(state == BACKGROUND);
                } else {
                    EXPECT_TRUE(state == BACKGROUND);
                }
            }
        }
    }
}


TEST(pore_id_tests, seal_gaps_cavity) {
    ProteinGrid g_cav = ProteinGrid("test_files/cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_cav);
    cylinder_completion(g_cav);
    perpendicularity_cylinder(g_cav);
    seal_gaps(g_cav);

    for (int x = g_cav.min.x; x < g_cav.max.x; x++) {
        for (int y = g_cav.min.y; y < g_cav.max.y; y++) {
            for (int z = g_cav.min.z; z < g_cav.max.z; z++) {
                uint8_t state = g_cav.at(x, y, z);
                if (state == OCCUPIED) {
                } else if (3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12) {
                    EXPECT_TRUE(state == POTENTIAL_PORE);
                } else {
                    EXPECT_TRUE(state == BACKGROUND);
                }
            }
        }
    }
}


TEST(pore_id_tests, identify_pore_nuclei) {
    ProteinGrid g_cav = ProteinGrid("test_files/cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_cav);
    cylinder_completion(g_cav);
    perpendicularity_cylinder(g_cav);
    seal_gaps(g_cav);
    remove_shallow_regions(g_cav, 1.0);
    identify_pore_nuclei(g_cav);

    for (int x = g_cav.min.x; x < g_cav.max.x; x++) {
        for (int y = g_cav.min.y; y < g_cav.max.y; y++) {
            for (int z = g_cav.min.z; z < g_cav.max.z; z++) {
                uint8_t state = g_cav.at(x, y, z);
                if (state == NUCLEUS) {
                    EXPECT_TRUE(x >= 7 && x <= 9 && y >= 7 && y <= 9 && z >= 7 && z <= 9);
                    EXPECT_TRUE(std::get<0>(g_cav.closest(x, y, z)) > g_cav.solvent_radius);
                }
            }
        }
    }
}

TEST(pore_id_tests, compute_volumes_cavity) {
    ProteinGrid g_cav = ProteinGrid("test_files/cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_cav);
    cylinder_completion(g_cav);
    perpendicularity_cylinder(g_cav);
    seal_gaps(g_cav);
    remove_shallow_regions(g_cav, 1.0);
    identify_pore_nuclei(g_cav);
    initial_pore_assembly(g_cav, 1.0);
    complete_pores(g_cav);
    label_pore_boxes(g_cav);
    lining_residues(g_cav);

    EXPECT_EQ(1, g_cav.clusters.size());
    EXPECT_EQ(1, g_cav.clusters[0].lining_residues.size());

    for (int x = g_cav.min.x; x < g_cav.max.x; x++) {
        for (int y = g_cav.min.y; y < g_cav.max.y; y++) {
            for (int z = g_cav.min.z; z < g_cav.max.z; z++) {
                uint8_t state = g_cav.at(x, y, z);
                if (3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12 && state != OCCUPIED) {
                    EXPECT_TRUE(state == IN_CLUSTER);
                    if (state != IN_CLUSTER) {
                        std::cout << Vec<int>(x, y, z) << std::endl;
                    }
                } else {
                    EXPECT_FALSE(state == IN_CLUSTER);
                }
            }
        }
    }
}

TEST(pore_id_tests, compute_volumes_pore) {
    ProteinGrid g_pore = ProteinGrid("test_files/pore.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_pore);
    cylinder_completion(g_pore);
    perpendicularity_cylinder(g_pore);
    seal_gaps(g_pore);
    remove_shallow_regions(g_pore, 1.0);
    identify_pore_nuclei(g_pore);
    initial_pore_assembly(g_pore, 1.0);
    complete_pores(g_pore);
    label_pore_boxes(g_pore);
    lining_residues(g_pore);

    EXPECT_EQ(1, g_pore.clusters.size());
    EXPECT_EQ(1, g_pore.clusters[0].lining_residues.size());

    for (int x = g_pore.min.x; x < g_pore.max.x; x++) {
        for (int y = g_pore.min.y; y < g_pore.max.y; y++) {
            for (int z = g_pore.min.z; z < g_pore.max.z; z++) {
                Vec<int> b = Vec<int>(x, y, z);
                uint8_t state = g_pore.at(x, y, z);

                if (state == OCCUPIED) {
                } else if (3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12) {
                    EXPECT_TRUE(state == IN_CLUSTER);
                } else if (b == Vec<int>(2, 6, 8) || b == Vec<int>(2, 7, 6) || b == Vec<int>(2, 7, 7)
                           || b == Vec<int>(2, 7, 8) || b == Vec<int>(2, 7, 9) || b == Vec<int>(2, 8, 6)
                           || b == Vec<int>(2, 8, 7) || b == Vec<int>(2, 8, 8) || b == Vec<int>(2, 8, 9)
                           || b == Vec<int>(2, 8, 10) || b == Vec<int>(2, 9, 7) || b == Vec<int>(2, 9, 8)
                           || b == Vec<int>(2, 9, 9) || b == Vec<int>(2, 10, 8) || b == Vec<int>(2, 10, 9)) {
                    EXPECT_FALSE(state == IN_CLUSTER);
                } else {
                    EXPECT_TRUE(state == BACKGROUND);
                }
            }
        }
    }
}

TEST(pore_id_tests, compute_volumes_pore_plus_cavity) {
    ProteinGrid g_pore = ProteinGrid("test_files/pore_plus_cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_pore);
    cylinder_completion(g_pore);
    perpendicularity_cylinder(g_pore);
    seal_gaps(g_pore);
    remove_shallow_regions(g_pore, 1.0);
    identify_pore_nuclei(g_pore);
    initial_pore_assembly(g_pore, 1.0);
    complete_pores(g_pore);
    label_pore_boxes(g_pore);
    lining_residues(g_pore);

    EXPECT_EQ(2, g_pore.clusters.size());
    EXPECT_EQ(1, g_pore.clusters[0].lining_residues.size());
    EXPECT_EQ(1, g_pore.clusters[1].lining_residues.size());
    EXPECT_TRUE(g_pore.clusters[0].pore);
    EXPECT_FALSE(g_pore.clusters[1].pore);

    for (int x = g_pore.min.x; x < g_pore.max.x; x++) {
        for (int y = g_pore.min.y; y < g_pore.max.y; y++) {
            for (int z = g_pore.min.z; z < g_pore.max.z; z++) {
                Vec<int> b = Vec<int>(x, y, z);
                BoxState state = g_pore.at(x, y, z);

                if (3 <= x && x <= 11 && 3 <= y && y <= 12 && 3 <= z && z <= 12 && state != OCCUPIED) {
                    if (b == Vec<int>(10, 6, 6)) {
                        EXPECT_TRUE(state < IN_CLUSTER);
                    } else {
                        EXPECT_TRUE(state >= IN_CLUSTER);
                        EXPECT_TRUE(std::find(g_pore.clusters[0].boxes.begin(), g_pore.clusters[0].boxes.end(),
                                              Vec<int>(x, y, z)) != g_pore.clusters[0].boxes.end());
                        EXPECT_FALSE(std::find(g_pore.clusters[1].boxes.begin(), g_pore.clusters[1].boxes.end(),
                                               Vec<int>(x, y, z)) != g_pore.clusters[1].boxes.end());

                    }
                } else if (12 <= x && x <= 21 && 3 <= y && y <= 12 && 3 <= z && z <= 12 && state != OCCUPIED) {
                    if (b == Vec<int>(12, 6, 8) || b == Vec<int>(12, 8, 6) || b == Vec<int>(12, 8, 10)
                        || b == Vec<int>(12, 10, 8)) {
                        EXPECT_TRUE(state < IN_CLUSTER);
                    } else {
                        EXPECT_TRUE(state >= IN_CLUSTER);
                        if (b == Vec<int>(12, 7, 8) || b == Vec<int>(12, 8, 7) || b == Vec<int>(12, 8, 8)
                            || b == Vec<int>(12, 8, 9) || b == Vec<int>(12, 9, 8)) {
                            EXPECT_FALSE(std::find(g_pore.clusters[1].boxes.begin(), g_pore.clusters[1].boxes.end(),
                                                   Vec<int>(x, y, z)) != g_pore.clusters[1].boxes.end());
                            EXPECT_TRUE(std::find(g_pore.clusters[0].boxes.begin(), g_pore.clusters[0].boxes.end(),
                                                  Vec<int>(x, y, z)) != g_pore.clusters[0].boxes.end());
                        } else {
                            EXPECT_TRUE(std::find(g_pore.clusters[1].boxes.begin(), g_pore.clusters[1].boxes.end(),
                                                  Vec<int>(x, y, z)) != g_pore.clusters[1].boxes.end());
                            EXPECT_FALSE(std::find(g_pore.clusters[0].boxes.begin(), g_pore.clusters[0].boxes.end(),
                                                   Vec<int>(x, y, z)) != g_pore.clusters[0].boxes.end());
                        }
                    }
                } else {
                    EXPECT_TRUE(state < IN_CLUSTER);
                }
            }
        }
    }
}

TEST(pore_id_tests, identify_pore_nuclei_standalone) {
    ProteinGrid g_cav = ProteinGrid("test_files/cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_cav);
    perpendicularity_standalone(g_cav);
    seal_gaps(g_cav);
    remove_shallow_regions(g_cav, 1.0);
    identify_pore_nuclei(g_cav);

    for (int x = g_cav.min.x; x < g_cav.max.x; x++) {
        for (int y = g_cav.min.y; y < g_cav.max.y; y++) {
            for (int z = g_cav.min.z; z < g_cav.max.z; z++) {
                uint8_t state = g_cav.at(x, y, z);
                if (state == NUCLEUS) {
                    EXPECT_TRUE(x >= 7 && x <= 9 && y >= 7 && y <= 9 && z >= 7 && z <= 9);
                    EXPECT_TRUE(std::get<0>(g_cav.closest(x, y, z)) > g_cav.solvent_radius);
                }
            }
        }
    }
}

TEST(pore_id_tests, compute_volumes_cavity_standalone) {
    ProteinGrid g_cav = ProteinGrid("test_files/cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_cav);
    perpendicularity_standalone(g_cav);
    seal_gaps(g_cav);
    remove_shallow_regions(g_cav, 1.0);
    identify_pore_nuclei(g_cav);
    initial_pore_assembly(g_cav, 1.0);
    complete_pores(g_cav);
    label_pore_boxes(g_cav);
    lining_residues(g_cav);

    EXPECT_EQ(1, g_cav.clusters.size());
    EXPECT_EQ(1, g_cav.clusters[0].lining_residues.size());

    for (int x = g_cav.min.x; x < g_cav.max.x; x++) {
        for (int y = g_cav.min.y; y < g_cav.max.y; y++) {
            for (int z = g_cav.min.z; z < g_cav.max.z; z++) {
                uint8_t state = g_cav.at(x, y, z);
                if (3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12 && state != OCCUPIED) {
                    EXPECT_TRUE(state == IN_CLUSTER);
                } else {
                    EXPECT_FALSE(state == IN_CLUSTER);
                }
            }
        }
    }
}

TEST(pore_id_tests, compute_volumes_pore_standalone) {
    ProteinGrid g_pore = ProteinGrid("test_files/pore.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_pore);
    perpendicularity_standalone(g_pore);
    seal_gaps(g_pore);
    remove_shallow_regions(g_pore, 1.0);
    identify_pore_nuclei(g_pore);
    initial_pore_assembly(g_pore, 1.0);
    complete_pores(g_pore);
    label_pore_boxes(g_pore);
    lining_residues(g_pore);

    EXPECT_EQ(1, g_pore.clusters.size());
    EXPECT_EQ(1, g_pore.clusters[0].lining_residues.size());
    EXPECT_TRUE(g_pore.clusters[0].pore);

    for (int x = g_pore.min.x; x < g_pore.max.x; x++) {
        for (int y = g_pore.min.y; y < g_pore.max.y; y++) {
            for (int z = g_pore.min.z; z < g_pore.max.z; z++) {
                Vec<int> b = Vec<int>(x, y, z);
                BoxState state = g_pore.at(x, y, z);

                if (state == OCCUPIED) {
                } else if (3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12) {
                    EXPECT_TRUE(state == IN_CLUSTER);
                } else {
                    EXPECT_TRUE(state < IN_CLUSTER);
                }
            }
        }
    }
}

TEST(pore_id_tests, compute_volumes_pore_plus_standalone) {
    ProteinGrid g_pore = ProteinGrid("test_files/pore_plus_cavity.pdb", "test_files/kept.pdb", "test_files/skipped.pdb", "test_files/log.yaml", 1.0, 1.0, false, false, true, true);
    collision_detection(g_pore);
    perpendicularity_standalone(g_pore);
    seal_gaps(g_pore);
    remove_shallow_regions(g_pore, 1.0);
    identify_pore_nuclei(g_pore);
    initial_pore_assembly(g_pore, 0.0);
    complete_pores(g_pore);
    label_pore_boxes(g_pore);
    lining_residues(g_pore);

    EXPECT_EQ(2, g_pore.clusters.size());
    EXPECT_EQ(1, g_pore.clusters[0].lining_residues.size());
    EXPECT_TRUE(g_pore.clusters[0].pore);
    EXPECT_EQ(1, g_pore.clusters[1].lining_residues.size());
    EXPECT_FALSE(g_pore.clusters[1].pore);
    EXPECT_TRUE(g_pore.clusters[0].pore);
    EXPECT_FALSE(g_pore.clusters[1].pore);
    EXPECT_TRUE(g_pore.clusters[0].size() >= 130);
    EXPECT_TRUE(g_pore.clusters[1].size() >= 110);

    for (int x = g_pore.min.x; x < g_pore.max.x; x++) {
        for (int y = g_pore.min.y; y < g_pore.max.y; y++) {
            for (int z = g_pore.min.z; z < g_pore.max.z; z++) {
                Vec<int> b = Vec<int>(x, y, z);
                BoxState state = g_pore.at(x, y, z);

                if (state >= IN_CLUSTER) {
                    EXPECT_TRUE(3 <= x && x <= 21 && 3 <= y && y <= 12 && 3 <= z && z <= 12);

                    if (std::find(g_pore.clusters[0].boxes.begin(), g_pore.clusters[0].boxes.end(),
                                  Vec<int>(x, y, z)) != g_pore.clusters[0].boxes.end()) {
                        EXPECT_TRUE(3 <= x && x <= 12 && 3 <= y && y <= 12 && 3 <= z && z <= 12);
                    } else if (std::find(g_pore.clusters[1].boxes.begin(), g_pore.clusters[1].boxes.end(),
                                         Vec<int>(x, y, z)) != g_pore.clusters[1].boxes.end()) {
                        EXPECT_TRUE(12 <= x && x <= 21 && 3 <= y && y <= 12 && 3 <= z && z <= 12);
                    }
                }
            }
        }
    }
}
