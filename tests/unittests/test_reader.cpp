#include <vector>
#include <memory>
#include "gtest/gtest.h"
#include "atom.h"
#include "enums.h"
#include "reader.h"


TEST(reader_tests, parse_PDB) {
    std::vector<std::shared_ptr<Atom>> atoms;
    parse_PDB("test_files/reader_test.pdb", atoms, false, false);
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
    EXPECT_EQ("EL", atoms[3]->element);
    EXPECT_EQ("+-", atoms[3]->charge);
    for (size_t i = 0; i < atoms.size(); i++) {
        EXPECT_EQ(i, atoms[i]->id);
    }
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", atoms, true, false);
    EXPECT_EQ(2, atoms.size());
    atoms.clear();
    parse_PDB("test_files/reader_test.pdb", atoms, false, true);
    EXPECT_EQ(5, atoms.size());
}

TEST(reader_tests, parse_PDB_simple_coord) {
    std::vector<std::shared_ptr<Atom>> atoms;
    parse_PDB("test_files/simple_coord.pdb", atoms, false, false);
    EXPECT_EQ(3, atoms.size());
}

TEST(reader_tests, pseudo_pdb) {
    EXPECT_EQ("ATOM  *****  C   AAA A   1     -13.267  87.421  34.691",
              pseudo_pdb(100000, Vec<double>(-13.267, 87.421, 34.691)));
    EXPECT_EQ("ATOM  99999  C   AAA A   1     -13.267  87.421  34.691",
              pseudo_pdb(99999, Vec<double>(-13.267, 87.421, 34.691)));
    EXPECT_EQ("ATOM      1  C   AAA A   1    -133.2688777.4211234.691",
              pseudo_pdb(1, Vec<double>(-133.2677, 8777.42122, 1234.69133)));
}

TEST(reader_tests, parse_and_write) {
    std::vector<std::shared_ptr<Atom>> atoms;
    parse_PDB("test_files/real.pdb", atoms, false, true);
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