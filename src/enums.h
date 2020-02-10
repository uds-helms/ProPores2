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
