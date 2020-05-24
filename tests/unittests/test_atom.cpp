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
#include "atom.h"
#include "enums.h"

/*
INPUT PDB
ATOM  *****  C   AAA A****     -13.267  87.421  34.691
ATOM      1  N   SER A   4     -13.267  87.421  34.691  1.00 50.42           N
ATOM      2  CA aGLU A   4i   -0.123451.1234562.1234560.123449.367           C +
ATOM      3  1H  SER A   4     -11.448  86.988  36.332  1.00 47.23          FE
*/


class AtomTest : public ::testing::Test {
public:
    AtomTest() :
            pseudo(std::make_shared<Atom>(Atom(0, "ATOM  *****  C   AAA A****     -13.267  87.421  34.691"))),
            element(std::make_shared<Atom>(
                    Atom(1, "ATOM      1  N   SER A   4     -13.267  87.421  34.691  1.00 50.42           N  "))),
            complete(std::make_shared<Atom>(
                    Atom(2, "ATOM      2  CA aGLU A   4i   -0.123451.1234562.1234560.123449.367           C +"))),
            simple(std::make_shared<Atom>(
                    Atom(3, "ATOM      3  1H  SER A   4     -11.448  86.988  36.332  1.00 47.23          FE  "))) {}

    std::shared_ptr<Atom> pseudo;
    std::shared_ptr<Atom> element;
    std::shared_ptr<Atom> complete;
    std::shared_ptr<Atom> simple;
};

TEST_F(AtomTest, constructor) {
    // pseudo
    EXPECT_EQ(0, pseudo->id);
    EXPECT_EQ(100000, pseudo->PDB_id);
    EXPECT_EQ(C, pseudo->type);
    EXPECT_EQ("C", pseudo->full_type);
    EXPECT_EQ("", pseudo->alternate_location);
    EXPECT_EQ(INVALID_RESIDUE, pseudo->residue_type);
    EXPECT_EQ(10000, pseudo->residue_id);
    EXPECT_EQ("A", pseudo->chain);
    EXPECT_EQ("", pseudo->insertion_code);
    EXPECT_EQ(Vec<double>(-13.267, 87.421, 34.691), pseudo->coord);
    EXPECT_EQ(Vec<double>(-13.267, 87.421, 34.691), pseudo->original);
    EXPECT_EQ("", pseudo->atom_occupancy);
    EXPECT_EQ("", pseudo->temperature);
    EXPECT_EQ(C, pseudo->element);
    EXPECT_EQ("", pseudo->charge);
    EXPECT_EQ(1.7, pseudo->vdw);
    EXPECT_EQ(12.01, pseudo->mass);
    EXPECT_EQ(nullptr, pseudo->KD_left);
    EXPECT_EQ(nullptr, pseudo->KD_right);
    EXPECT_EQ(false, pseudo->in_tree);
    // with element at the end
    EXPECT_EQ(1, element->id);
    EXPECT_EQ(1, element->PDB_id);
    EXPECT_EQ(N, element->type);
    EXPECT_EQ("N", element->full_type);
    EXPECT_EQ("", element->alternate_location);
    EXPECT_EQ(SER, element->residue_type);
    EXPECT_EQ(4, element->residue_id);
    EXPECT_EQ("A", element->chain);
    EXPECT_EQ("", element->insertion_code);
    EXPECT_EQ(Vec<double>(-13.267, 87.421, 34.691), element->coord);
    EXPECT_EQ(Vec<double>(-13.267, 87.421, 34.691), element->original);
    EXPECT_EQ("1.00", element->atom_occupancy);
    EXPECT_EQ("50.42", element->temperature);
    EXPECT_EQ(N, element->element);
    EXPECT_EQ("", element->charge);
    EXPECT_EQ(1.55, element->vdw);
    EXPECT_EQ(14.01, element->mass);
    EXPECT_EQ(nullptr, element->KD_left);
    EXPECT_EQ(nullptr, element->KD_right);
    EXPECT_EQ(false, element->in_tree);
    // complete information
    EXPECT_EQ(2, complete->id);
    EXPECT_EQ(2, complete->PDB_id);
    EXPECT_EQ(C, complete->type);
    EXPECT_EQ("CA", complete->full_type);
    EXPECT_EQ("a", complete->alternate_location);
    EXPECT_EQ(GLU, complete->residue_type);
    EXPECT_EQ(4, complete->residue_id);
    EXPECT_EQ("A", complete->chain);
    EXPECT_EQ("i", complete->insertion_code);
    EXPECT_EQ(Vec<double>(-0.12345, 1.123456, 2.123456), complete->coord);
    EXPECT_EQ(Vec<double>(-0.12345, 1.123456, 2.123456), complete->original);
    EXPECT_EQ("0.1234", complete->atom_occupancy);
    EXPECT_EQ("49.367", complete->temperature);
    EXPECT_EQ(C, complete->element);
    EXPECT_EQ("+", complete->charge);
    EXPECT_EQ(1.7, complete->vdw);
    EXPECT_EQ(12.01, complete->mass);
    EXPECT_EQ(nullptr, complete->KD_left);
    EXPECT_EQ(nullptr, complete->KD_right);
    EXPECT_EQ(false, complete->in_tree);
    // simple
    EXPECT_EQ(3, simple->id);
    EXPECT_EQ(3, simple->PDB_id);
    EXPECT_EQ(Fe, simple->type);
    EXPECT_EQ("", simple->alternate_location);
    EXPECT_EQ(SER, simple->residue_type);
    EXPECT_EQ(4, simple->residue_id);
    EXPECT_EQ("A", simple->chain);
    EXPECT_EQ("", simple->insertion_code);
    EXPECT_EQ(Vec<double>(-11.448, 86.988, 36.332), simple->coord);
    EXPECT_EQ(Vec<double>(-11.448, 86.988, 36.332), simple->original);
    EXPECT_EQ("1.00", simple->atom_occupancy);
    EXPECT_EQ("47.23", simple->temperature);
    EXPECT_EQ(Fe, simple->element);
    EXPECT_EQ("", simple->charge);
    EXPECT_EQ(1.96, simple->vdw);
    EXPECT_EQ(55.85, simple->mass);
    EXPECT_EQ(nullptr, simple->KD_left);
    EXPECT_EQ(nullptr, simple->KD_right);
    EXPECT_EQ(false, simple->in_tree);
}

TEST_F(AtomTest, KD) {
    pseudo->KD_left = complete;
    EXPECT_EQ(true, pseudo->has_left());
    EXPECT_EQ(false, pseudo->has_right());
    EXPECT_EQ(2, pseudo->KD_left->id);
    pseudo->KD_right = simple;
    EXPECT_EQ(true, pseudo->has_left());
    EXPECT_EQ(true, pseudo->has_right());
    EXPECT_EQ(3, pseudo->KD_right->id);
}

TEST(atom_tests, vdw_radius) {
    EXPECT_EQ(1.55, vdw_radius(N));
    EXPECT_EQ(1.7, vdw_radius(C));
    EXPECT_EQ(1.52, vdw_radius(O));
    EXPECT_EQ(1.2, vdw_radius(H));
    EXPECT_EQ(1.8, vdw_radius(S));
    EXPECT_EQ(0.0, vdw_radius(INVALID_ATOM));
}

TEST(atom_tests, atomic_mass) {
    EXPECT_EQ(14.01, atomic_mass(N));
    EXPECT_EQ(12.01, atomic_mass(C));
    EXPECT_EQ(16.00, atomic_mass(O));
    EXPECT_EQ(1.008, atomic_mass(H));
    EXPECT_EQ(32.07, atomic_mass(S));
    EXPECT_EQ(0.0, atomic_mass(INVALID_ATOM));
}

TEST(atom_tests, parse_atom_type) {
    std::pair<AtomType, std::string> res = std::make_pair(H, "1HB");
    EXPECT_EQ(res, parse_atom_type("1HB ", ATOM));
    EXPECT_EQ(res, parse_atom_type(" 1HB", ATOM));
    EXPECT_EQ(res, parse_atom_type("1HB ", HETATOM));
    res = std::make_pair(H, "21HB");
    EXPECT_EQ(res, parse_atom_type("21HB ", ATOM));
    EXPECT_EQ(res, parse_atom_type("21HB ", HETATOM));
    res = std::make_pair(C, "CA");
    EXPECT_EQ(res, parse_atom_type("CA  ", ATOM));
    EXPECT_EQ(res, parse_atom_type(" CA ", ATOM));
    EXPECT_EQ(res, parse_atom_type("  CA", ATOM));
    EXPECT_EQ(res, parse_atom_type(" CA ", HETATOM));
    res = std::make_pair(Ca, "CA");
    EXPECT_EQ(res, parse_atom_type("CA  ", HETATOM));
    EXPECT_EQ(res, parse_atom_type("  CA", HETATOM));
    res = std::make_pair(I, "I");
    EXPECT_EQ(res, parse_atom_type("I   ", ATOM));
    EXPECT_EQ(res, parse_atom_type("   I", ATOM));
    EXPECT_EQ(res, parse_atom_type(" I  ", ATOM));
    EXPECT_EQ(res, parse_atom_type("  I ", ATOM));
    EXPECT_EQ(res, parse_atom_type("I   ", HETATOM));
    EXPECT_EQ(res, parse_atom_type("   I", HETATOM));
    EXPECT_EQ(res, parse_atom_type(" I  ", HETATOM));
    EXPECT_EQ(res, parse_atom_type("  I ", HETATOM));
    res = std::make_pair(INVALID_ATOM, "");
    EXPECT_EQ(res, parse_atom_type("1C", ATOM));
    EXPECT_EQ(res, parse_atom_type("1CH", ATOM));
    EXPECT_EQ(res, parse_atom_type("12HBB", ATOM));
    EXPECT_EQ(res, parse_atom_type("", ATOM));
}

TEST_F(AtomTest, sort_atoms) {
    std::vector<std::shared_ptr<Atom>> atoms = {complete, pseudo, simple, element};
    sort_atoms(atoms);
    EXPECT_EQ(0, atoms[0]->id);
    EXPECT_EQ(1, atoms[1]->id);
    EXPECT_EQ(2, atoms[2]->id);
    EXPECT_EQ(3, atoms[3]->id);
}
