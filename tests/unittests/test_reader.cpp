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

#include <vector>
#include <memory>
#include "gtest/gtest.h"
#include "atom.h"
#include "enums.h"
#include "reader.h"


TEST(reader_tests, parse_PDB) {
    std::vector<std::shared_ptr<Atom>> atoms;
    PDBReaderStats stats = PDBReaderStats();
    // all
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, REMOVE_ALL_HETERO_ATOMS, false, true);
    EXPECT_EQ(4, atoms.size());
    EXPECT_EQ(12345, atoms[3]->PDB_id);
    EXPECT_EQ(C, atoms[3]->type);
    EXPECT_EQ("A", atoms[3]->alternate_location);
    EXPECT_EQ(ALA, atoms[3]->residue_type);
    EXPECT_EQ("A", atoms[3]->chain);
    EXPECT_EQ(1234, atoms[3]->residue_id);
    EXPECT_EQ("i", atoms[3]->insertion_code);
    EXPECT_EQ("0.0001", atoms[3]->atom_occupancy);
    EXPECT_EQ("123456", atoms[3]->temperature);
    EXPECT_EQ(C, atoms[3]->element);
    EXPECT_EQ("+-", atoms[3]->charge);
    for (size_t i = 0; i < atoms.size(); i++) {
        EXPECT_EQ(i, atoms[i]->id);
    }
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, KEEP_ALL_HETERO_ATOMS, false, true);
    EXPECT_EQ(7, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, REMOVE_ALL_H_ATOMS, REMOVE_ALL_HETERO_ATOMS, false, true);

    EXPECT_EQ(2, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, KEEP_ALL_HETERO_ATOMS, true, false);
    EXPECT_EQ(13, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, REMOVE_ONLY_PROTEIN_H_ATOMS, KEEP_ALL_HETERO_ATOMS, true, false);
    EXPECT_EQ(10, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, REMOVE_ONLY_HETERO_H_ATOMS, KEEP_ALL_HETERO_ATOMS, true, false);
    EXPECT_EQ(12, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, REMOVE_ALL_H_ATOMS, KEEP_ALL_HETERO_ATOMS, true, false);
    EXPECT_EQ(9, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, REMOVE_ALL_HETERO_ATOMS, true, false);
    EXPECT_EQ(10, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, REMOVE_HETERO_ATOMS_EXCEPT_DUMMY, true, false);
    EXPECT_EQ(11, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, REMOVE_ONLY_DUMMY_HETERO_ATOMS, true, false);
    EXPECT_EQ(12, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, REMOVE_ONLY_DUMMY_HETERO_ATOMS, true, true);
    EXPECT_EQ(7, atoms.size());
}

TEST(reader_tests, parse_PDB_simple_coord) {
    std::vector<std::shared_ptr<Atom>> atoms;
    PDBReaderStats stats = PDBReaderStats();
    parse_PDB("test_files/simple_coord.pdb", "test_files/reader_output.pdb", "test_files/reader_output.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, KEEP_ALL_HETERO_ATOMS, false, false);
    EXPECT_EQ(3, atoms.size());
}

TEST(reader_tests, pseudo_pdb) {
    EXPECT_EQ("ATOM  *****  C   AAA A   1     -13.267  87.421  34.691  1.00  0.00           C  ",
              pseudo_pdb(100000, Vec<double>(-13.267, 87.421, 34.691)));
    EXPECT_EQ("ATOM  99999  C   AAA A   1     -13.267  87.421  34.691  1.00  0.00           C  ",
              pseudo_pdb(99999, Vec<double>(-13.267, 87.421, 34.691)));
    EXPECT_EQ("ATOM      1  C   AAA A   1    -133.2688777.4211234.691  1.00  0.00           C  ",
              pseudo_pdb(1, Vec<double>(-133.2677, 8777.42122, 1234.69133)));
}

TEST(reader_tests, parse_and_write) {
    std::vector<std::shared_ptr<Atom>> atoms;
    PDBReaderStats stats = PDBReaderStats();
    parse_PDB("test_files/real.pdb", "test_files/kept.pdb", "test_files/skipped.pdb",
              atoms, stats, KEEP_ALL_H_ATOMS, KEEP_ALL_HETERO_ATOMS, true, false);
    EXPECT_EQ(4389, atoms.size());

    std::ifstream real_file("test_files/real.pdb");
    std::vector<std::string> real_lines;
    std::string line;
    while (std::getline(real_file, line)) { real_lines.push_back(r_strip(line)); }
    real_file.close();

    write_PDB("test_files/real_test_output.pdb", atoms);
    std::ifstream test_file("test_files/real_test_output.pdb");
    std::vector<std::string> test_lines;
    while (std::getline(test_file, line)) { test_lines.push_back(r_strip(line)); }
    test_file.close();

    EXPECT_EQ(real_lines.size(), test_lines.size());
    for (size_t i = 0; i < real_lines.size(); i++) {
        EXPECT_EQ(real_lines[i], test_lines[i]);
    }
}
