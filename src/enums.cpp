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

#include <string>
#include <vector>
#include "enums.h"

std::string to_str(const AtomType atom) {
    switch (atom) {
        case H:
            return "H";
        case O:
            return "O";
        case N:
            return "N";
        case S:
            return "S";
        case C:
            return "C";
        default:
            return "INVALID_ATOM";
    }
}

AtomType to_atom_type(const std::string &str) {
    if (str == "H") return H;
    if (str == "O") return O;
    if (str == "N") return N;
    if (str == "S") return S;
    if (str == "C") return C;
    return INVALID_ATOM;
}


std::string to_str(const ResidueType residue) {
    switch (residue) {
        case ALA:
            return "ALA";
        case ARG:
            return "ARG";
        case ASN:
            return "ASN";
        case ASP:
            return "ASP";
        case CYS:
            return "CYS";
        case GLN:
            return "GLN";
        case GLU:
            return "GLU";
        case GLY:
            return "GLY";
        case HIS:
            return "HIS";
        case ILE:
            return "ILE";
        case LEU:
            return "LEU";
        case LYS:
            return "LYS";
        case MET:
            return "MET";
        case PHE:
            return "PHE";
        case PRO:
            return "PRO";
        case SER:
            return "SER";
        case THR:
            return "THR";
        case TRP:
            return "TRP";
        case TYR:
            return "TYR";
        case VAL:
            return "VAL";
        default:
            return "INVALID_RESIDUE";
    }
}

ResidueType to_residue_type(const std::string &str) {
    if (str == "ALA") return ALA;
    if (str == "ARG") return ARG;
    if (str == "ASN") return ASN;
    if (str == "ASP") return ASP;
    if (str == "CYS") return CYS;
    if (str == "GLN") return GLN;
    if (str == "GLU") return GLU;
    if (str == "GLY") return GLY;
    if (str == "HIS") return HIS;
    if (str == "ILE") return ILE;
    if (str == "LEU") return LEU;
    if (str == "LYS") return LYS;
    if (str == "MET") return MET;
    if (str == "PHE") return PHE;
    if (str == "PRO") return PRO;
    if (str == "SER") return SER;
    if (str == "THR") return THR;
    if (str == "TRP") return TRP;
    if (str == "TYR") return TYR;
    if (str == "VAL") return VAL;
    return INVALID_RESIDUE;
}

std::string to_str(const BoxState state) {
    switch (state) {
        case OCCUPIED:
            return "OCCUPIED";
        case BACKGROUND:
            return "BACKGROUND";
        case POTENTIAL_PORE:
            return "POTENTIAL_PORE";
        case NUCLEUS:
            return "NUCLEUS";
        case IN_CLUSTER:
            return "IN_CLUSTER";
        case ON_CLUSTER_INSIDE:
            return "ON_CLUSTER_INSIDE";
        case TOUCHES_BACKGROUND:
            return "TOUCHES_BACKGROUND";
        case TOUCHES_PROTEIN:
            return "TOUCHES_PROTEIN";
        case TOUCHES_OTHER_CLUSTER:
            return "TOUCHES_OTHER_CLUSTER";
        case TOUCHES_POTENTIAL_PORE:
            return "TOUCHES_POTENTIAL_PORE";
        default:
            return "UNCLASSIFIED";
    }
}

BoxState to_box_state(const std::string &str) {
    if (str == "OCCUPIED") return OCCUPIED;
    if (str == "BACKGROUND") return BACKGROUND;
    if (str == "POTENTIAL_PORE") return POTENTIAL_PORE;
    if (str == "NUCLEUS") return NUCLEUS;
    if (str == "IN_CLUSTER") return IN_CLUSTER;
    if (str == "ON_CLUSTER_INSIDE") return ON_CLUSTER_INSIDE;
    if (str == "TOUCHES_BACKGROUND") return TOUCHES_BACKGROUND;
    if (str == "TOUCHES_PROTEIN") return TOUCHES_PROTEIN;
    if (str == "TOUCHES_OTHER_CLUSTER") return TOUCHES_OTHER_CLUSTER;
    if (str == "TOUCHES_POTENTIAL_PORE") return TOUCHES_POTENTIAL_PORE;
    return UNCLASSIFIED;
}

std::string to_str(const GateDifficulty difficulty) {
    switch (difficulty) {
        case EASY:
            return "EASY";
        case MEDIUM:
            return "MEDIUM";
        case HARD:
            return "HARD";
        default:
            return "DIFFICULTY_ERROR";
    }
}

GateDifficulty to_gate_difficulty(const std::string &str) {
    if (str == "EASY") return EASY;
    if (str == "MEDIUM") return MEDIUM;
    if (str == "HARD") return HARD;
    return DIFFICULTY_ERROR;
}

// map residue types to backbone definitions
std::vector<std::string> residue_definition(const ResidueType type) {
    if (type == ARG) return {"N", "CA", "CB", "CG", "CD", "NE", "CZ"};
    if (type == ASN) return {"N", "CA", "CB", "CG", "OD1"};
    if (type == ASP) return {"N", "CA", "CB", "CG", "OD1"};
    if (type == CYS) return {"N", "CA", "CB", "SG"};
    if (type == GLN) return {"N", "CA", "CB", "CG", "CD", "OE1"};
    if (type == GLU) return {"N", "CA", "CB", "CG", "CD", "OE1"};
    if (type == HIS) return {"N", "CA", "CB", "CG", "ND1"};
    if (type == ILE) return {"N", "CA", "CB", "CG1", "CD1"};
    if (type == LEU) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == LYS) return {"N", "CA", "CB", "CG", "CD", "CE", "NZ"};
    if (type == MET) return {"N", "CA", "CB", "CG", "SD", "CE"};
    if (type == PHE) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == PRO) return {"N", "CA", "CB", "CG", "CD"};
    if (type == SER) return {"N", "CA", "CB", "OG"};
    if (type == THR) return {"N", "CA", "CB", "OG1"};
    if (type == TRP) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == TYR) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == VAL) return {"N", "CA", "CB", "CG1"};
    return {};
}