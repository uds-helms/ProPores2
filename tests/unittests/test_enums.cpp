#include <vector>
#include <string>
#include "gtest/gtest.h"
#include "enums.h"

TEST(enum_tests, to_atom_type) {
    for (const std::string &str: {"H", "C", "O", "N", "S"}) {
        EXPECT_FALSE(INVALID_ATOM == to_atom_type(str));
    }
    for (const std::string &str: {"", "E", "  ", "1", "H1", "CA", "h"}) {
        EXPECT_TRUE(INVALID_ATOM == to_atom_type(str));
    }
}

TEST(enum_tests, atom_type_to_str) {
    std::vector<std::string> expected = {"H", "O", "N", "S", "C", "INVALID_ATOM"};
    std::vector<std::string> actual;
    for (const AtomType &type: {H, O, N, S, C, INVALID_ATOM}) { actual.push_back(to_str(type)); }
    EXPECT_EQ(expected, actual);
}

TEST(enum_tests, to_residue_type) {
    for (const std::string &str: {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                                  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}) {
        EXPECT_FALSE(INVALID_RESIDUE == to_residue_type(str));
    }
    for (const std::string &str: {"", "AL", "  ", "1", "ala", "blablub"}) {
        EXPECT_TRUE(INVALID_RESIDUE == to_residue_type(str));
    }
}

TEST(enum_tests, residue_type_to_str) {
    std::vector<std::string> expected = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                                         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                                         "INVALID_RESIDUE"};
    std::vector<std::string> actual;
    for (const ResidueType &type: {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR,
                                   TRP, TYR, VAL, INVALID_RESIDUE}) { actual.push_back(to_str(type)); }
    EXPECT_EQ(expected, actual);
}

TEST(enum_tests, to_box_state) {
    for (const std::string &str: {"OCCUPIED", "BACKGROUND", "POTENTIAL_PORE", "NUCLEUS", "IN_CLUSTER",
                                  "ON_CLUSTER_INSIDE", "TOUCHES_BACKGROUND", "TOUCHES_PROTEIN", "TOUCHES_OTHER_CLUSTER",
                                  "TOUCHES_POTENTIAL_PORE"}) {
        EXPECT_FALSE(UNCLASSIFIED == to_box_state(str));
    }
    for (const std::string &str: {"", "AL", "  ", "1", "ala", "blablub", "background"}) {
        EXPECT_TRUE(UNCLASSIFIED == to_box_state(str));
    }
}

TEST(enum_tests, box_state_to_str) {
    std::vector<std::string> expected = {"OCCUPIED", "BACKGROUND", "POTENTIAL_PORE", "NUCLEUS", "IN_CLUSTER",
                                         "ON_CLUSTER_INSIDE", "TOUCHES_BACKGROUND", "TOUCHES_PROTEIN",
                                         "TOUCHES_OTHER_CLUSTER", "TOUCHES_POTENTIAL_PORE", "UNCLASSIFIED"};
    std::vector<std::string> actual;
    for (const BoxState &diff: {OCCUPIED, BACKGROUND, POTENTIAL_PORE, NUCLEUS, IN_CLUSTER, ON_CLUSTER_INSIDE,
                                TOUCHES_BACKGROUND, TOUCHES_PROTEIN, TOUCHES_OTHER_CLUSTER,
                                TOUCHES_POTENTIAL_PORE, UNCLASSIFIED}) { actual.push_back(to_str(diff)); }
    EXPECT_EQ(expected, actual);
}

TEST(enum_tests, to_difficulty) {
    for (const std::string &str: {"EASY", "MEDIUM", "HARD"}) {
        EXPECT_FALSE(DIFFICULTY_ERROR == to_gate_difficulty(str));
    }
    for (const std::string &str: {"", "AL", "  ", "1", "ala", "blablub", "easy", "medium", "hard"}) {
        EXPECT_TRUE(DIFFICULTY_ERROR == to_gate_difficulty(str));
    }
}

TEST(enum_tests, difficulty_to_str) {
    std::vector<std::string> expected = {"EASY", "MEDIUM", "HARD", "DIFFICULTY_ERROR"};
    std::vector<std::string> actual;
    for (const GateDifficulty &diff: {EASY, MEDIUM, HARD, DIFFICULTY_ERROR}) { actual.push_back(to_str(diff)); }
    EXPECT_EQ(expected, actual);
}

