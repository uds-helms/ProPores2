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

#ifndef PROPORES_ENUMS_H
#define PROPORES_ENUMS_H

#include <string>

enum AtomType : uint8_t {
    INVALID_ATOM,
    H,
    O,
    N,
    S,
    C
};

enum ResidueType : uint8_t {
    INVALID_RESIDUE,
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLN,
    GLU,
    GLY,
    HIS,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL
};

enum BoxState : uint8_t {
    UNCLASSIFIED,
    OCCUPIED,
    BACKGROUND,
    POTENTIAL_PORE,
    NUCLEUS,
    IN_CLUSTER,
    ON_CLUSTER_INSIDE,
    TOUCHES_BACKGROUND,
    TOUCHES_PROTEIN,
    TOUCHES_OTHER_CLUSTER,
    TOUCHES_POTENTIAL_PORE
};

enum RemoveTag : uint8_t {
    REMOVE_NOTHING,
    REMOVE_CAVITIES,
    REMOVE_PORES
};

enum CylinderTag : uint8_t {
    AUTODETECT,
    RAY_TRACE,
    STANDALONE
};

enum GateDifficulty : uint8_t {
    EASY,
    MEDIUM,
    HARD,
    DIFFICULTY_ERROR
};


// convert enums to string
std::string to_str(AtomType atom);

std::string to_str(ResidueType residue);

std::string to_str(BoxState state);

std::string to_str(GateDifficulty difficulty);

// convert strings to enums
AtomType to_atom_type(const std::string &str);

ResidueType to_residue_type(const std::string &str);

BoxState to_box_state(const std::string &str);

GateDifficulty to_gate_difficulty(const std::string &str);

// map residue types to backbone definitions
std::vector<std::string> residue_definition(ResidueType type);

#endif //PROPORES_ENUMS_H
