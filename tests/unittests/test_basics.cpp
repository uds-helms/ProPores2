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
#include "basics.h"
#include <vector>

TEST(basics_tests, remove_spaces) {
    EXPECT_EQ("", remove_spaces("    "));
    EXPECT_EQ("", remove_spaces(""));
    EXPECT_EQ("1", remove_spaces(" 1"));
    EXPECT_EQ("1", remove_spaces("1 "));
    EXPECT_EQ("1", remove_spaces(" 1 "));
    EXPECT_EQ("1", remove_spaces("1"));
    EXPECT_EQ("blablub", remove_spaces("\t bla blub \n"));
}

TEST(basics_tests, relative_search_box) {
    RelativeSearchBox box_1 = RelativeSearchBox(0, 0, 0);
    EXPECT_EQ(Vec<int>(0, 0, 0), box_1.grid);
    EXPECT_EQ(0.0, box_1.length_sq);
    RelativeSearchBox box_2 = RelativeSearchBox(1, 2, 3);
    EXPECT_EQ(Vec<int>(1, 2, 3), box_2.grid);
    EXPECT_EQ(14.0, box_2.length_sq);
}

TEST(basics_tests, relative_search_spaces) {
    RelativeSearchSpaces spaces = RelativeSearchSpaces();
    EXPECT_TRUE(spaces.space_map.empty());
    EXPECT_EQ(1, spaces[0].size());
    EXPECT_EQ(pow(3, 3), spaces[1].size());
    EXPECT_EQ(pow(5, 3), spaces[2].size());
    EXPECT_EQ(3, spaces.space_map.size());
    EXPECT_EQ(Vec<int>(-1, -1, -1), spaces[1][0].grid);
    EXPECT_EQ(Vec<int>(1, 1, 1), spaces[1][pow(3, 3) - 1].grid);
}

TEST(basics_tests, split) {
    std::string str = "example\tstring";
    std::vector<std::string> res = {"example", "string"};
    EXPECT_EQ(res, split(str, '\t'));
    EXPECT_EQ(res, split(str));
    str = " example \t string ";
    EXPECT_EQ(res, split(str));
    res = {" example ", " string "};
    EXPECT_EQ(res, split(str, '\t'));
    str = "example";
    res = {"example"};
    EXPECT_EQ(res, split(str, '\t'));
    EXPECT_EQ(res, split(str));
    str = " example ";
    EXPECT_EQ(res, split(str));
    res = {" example "};
    EXPECT_EQ(res, split(str, '\t'));
    str = "";
    res = {};
    EXPECT_EQ(res, split(str, '\t'));
    EXPECT_EQ(res, split(str));
}

TEST(basics_tests, indent) {
    EXPECT_EQ("", indent(0));
    EXPECT_EQ("  ", indent(1));
    EXPECT_EQ("    ", indent(2));
}

TEST(basics_tests, r_strip) {
    EXPECT_EQ("  example", r_strip("  example \t\n "));
    EXPECT_EQ("", r_strip("\n   \t  "));
    EXPECT_EQ("bla blub", r_strip("bla blub"));
    EXPECT_EQ("", r_strip(""));
}

TEST(basics_tests, substr) {
    std::string str = "123456789";
    EXPECT_EQ("", substr("", 0, 0));
    EXPECT_EQ("", substr("", 0, 0));
    EXPECT_EQ("", substr(str, 0, 0));
    EXPECT_EQ("1", substr(str, 0, 1));
    EXPECT_EQ("12", substr(str, 0, 2));
    EXPECT_EQ("9", substr(str, 8, 5));
    EXPECT_EQ("", substr(str, 9, 5));
}

TEST(basics_tests, pdb_substr) {
    std::string str = "123456789***";
    EXPECT_EQ("", pdb_substr("", 0, 0));
    EXPECT_EQ("", pdb_substr("", 0, 0));
    EXPECT_EQ("", pdb_substr(str, 0, 0));
    EXPECT_EQ("1", pdb_substr(str, 0, 1));
    EXPECT_EQ("12", pdb_substr(str, 0, 2));
    EXPECT_EQ("9", pdb_substr(str, 8, 5));
    EXPECT_EQ("", pdb_substr(str, 9, 5));
    EXPECT_EQ("", pdb_substr(" ** *", 0, 10));
    EXPECT_EQ(100.0, pdb_substr(" ** *", 0, 10, 100.0));
    EXPECT_EQ(-2.3, pdb_substr(" ** *", 0, 10, -2.3));
}

TEST(basics_tests, wrap) {
    std::string str = "  12345    123456789   123 \t 12\n123 123";
    std::vector<std::string> expected = {"12345", "123456789", "123 12", "123", "123"};
    EXPECT_EQ(expected, wrap(str, 6));
}