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

enum RecordType : uint8_t {
    INVALID_RECORD,
    ATOM_RECORD,
    HETERO_RECORD
};

enum AtomType : uint8_t {
    INVALID_ATOM,
    // protein atoms
    C,
    H,
    N,
    O,
    S,
    // hetero atoms
    Ac,
    Ag,
    Al,
    Am,
    Ar,
    As,
    At,
    Au,
    B,
    Ba,
    Be,
    Bh,
    Bi,
    Bk,
    Br,
    Ca,
    Cd,
    Ce,
    Cf,
    Cl,
    Cm,
    Cn,
    Co,
    Cr,
    Cs,
    Cu,
    Db,
    Ds,
    Dy,
    Er,
    Es,
    Eu,
    F,
    Fe,
    Fl,
    Fm,
    Fr,
    Ga,
    Gd,
    Ge,
    He,
    Hf,
    Hg,
    Ho,
    Hs,
    I,
    In,
    Ir,
    K,
    Kr,
    La,
    Li,
    Lr,
    Lu,
    Lv,
    Mc,
    Md,
    Mg,
    Mn,
    Mo,
    Mt,
    Na,
    Nb,
    Nd,
    Ne,
    Nh,
    Ni,
    No,
    Np,
    Og,
    Os,
    P,
    Pa,
    Pb,
    Pd,
    Pm,
    Po,
    Pr,
    Pt,
    Pu,
    Ra,
    Rb,
    Re,
    Rf,
    Rg,
    Rh,
    Rn,
    Ru,
    Sb,
    Sc,
    Se,
    Sg,
    Si,
    Sm,
    Sn,
    Sr,
    Ta,
    Tb,
    Tc,
    Te,
    Th,
    Ti,
    Tl,
    Tm,
    Ts,
    U,
    V,
    W,
    Xe,
    Y,
    Yb,
    Zn,
    Zr
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

enum HAtomTag : uint8_t {
    KEEP_ALL_H_ATOMS,
    REMOVE_ALL_H_ATOMS,
    REMOVE_ONLY_PROTEIN_H_ATOMS,
    REMOVE_ONLY_HETERO_H_ATOMS,
    INVALID_H_ATOM_TAG
};

enum HeteroAtomTag : uint8_t {
    KEEP_ALL_HETERO_ATOMS,
    REMOVE_ALL_HETERO_ATOMS,
    REMOVE_HETERO_ATOMS_EXCEPT_DUMMY,
    REMOVE_ONLY_DUMMY_HETERO_ATOMS,
    INVALID_HETERO_ATOM_TAG
};

enum RemoveTag : uint8_t {
    REMOVE_NOTHING,
    REMOVE_CAVITIES,
    REMOVE_PORES,
    INVALID_REMOVE_TAG
};

enum ComputationMode : uint8_t {
    AUTODETECT,
    RAY_TRACE,
    STANDALONE,
    INVALID_COMPUTATION_MODE
};

enum GateDifficulty : uint8_t {
    EASY,
    MEDIUM,
    HARD,
    DIFFICULTY_ERROR
};

enum PreparationTag: uint8_t {
    AXIS_AND_GATE,
    ONLY_AXIS,
    ONLY_GATE,
    NO_PREPARATION,
    INVALID_PREPARATION_TAG
};


// convert enums to string
std::string to_str(RecordType record);
std::string to_str(AtomType atom);
std::string to_str(ResidueType residue);
std::string to_str(BoxState state);
std::string to_str(HAtomTag tag);
std::string to_str(HeteroAtomTag tag);
std::string to_str(GateDifficulty difficulty);
std::string to_str(ComputationMode tag);
std::string to_str(RemoveTag tag);
std::string to_str(PreparationTag tag);

// convert strings to enums
RecordType to_record_type(const std::string &str);
AtomType to_atom_type(const std::string &str);
AtomType to_atom_type(const char &str);
ResidueType to_residue_type(const std::string &str);
BoxState to_box_state(const std::string &str);
HAtomTag to_h_atom_tag(const std::string &str);
HeteroAtomTag to_hetero_atom_tag(const std::string &str);
GateDifficulty to_gate_difficulty(const std::string &str);
ComputationMode to_computation_mode(const std::string &str);
RemoveTag to_remove_tag(const std::string &str);
PreparationTag to_preparation_tag(const std::string &str);

// map residue types to backbone definitions
std::vector<std::string> residue_definition(ResidueType type);

#endif //PROPORES_ENUMS_H
